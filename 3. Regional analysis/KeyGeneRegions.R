library(vcfR)
library(tidyverse)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: Load metadata
# ══════════════════════════════════════════════════════════════════════════════

metadata <- read.csv("Heteroplasmy_analysis/Merged_cohorts/sample_data.csv", header=TRUE) %>%
  dplyr::select(sample_id, ancestry, macrohaplogroup)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Extract GT and AF from vcf_filtered
# ══════════════════════════════════════════════════════════════════════════════

gt_matrix     <- extract.gt(vcf_filtered, element = "GT", as.numeric = FALSE)
af_matrix     <- extract.gt(vcf_filtered, element = "AF", as.numeric = FALSE)
vcf_positions <- as.integer(vcf_filtered@fix[, "POS"])
vcf_ref       <- vcf_filtered@fix[, "REF"]
vcf_alt       <- vcf_filtered@fix[, "ALT"]
vcf_samples   <- colnames(gt_matrix)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: Extract heteroplasmic calls into long-format data frame
# ══════════════════════════════════════════════════════════════════════════════

het_logical <- gt_matrix == "0/1" | gt_matrix == "0/2"
het_calls   <- which(het_logical, arr.ind = TRUE)
cat("Total heteroplasmic calls:", nrow(het_calls), "\n")

# Clean AF values
clean_af <- function(x) {
  if (is.na(x) || x == "." || x == "") return(NA)
  as.numeric(x)
}

het_df <- data.frame(
  sample_id = vcf_samples[het_calls[, 2]],
  Position  = vcf_positions[het_calls[, 1]],
  REF       = vcf_ref[het_calls[, 1]],
  ALT       = vcf_alt[het_calls[, 1]],
  GT        = gt_matrix[het_calls],
  AF        = sapply(af_matrix[het_calls], clean_af)
) %>%
  left_join(metadata, by = "sample_id")

cat("Samples with heteroplasmic calls:", n_distinct(het_df$sample_id), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: Assign mitochondrial regions
#   12S rRNA : 577  – 1601
#   ND5      : 12337 – 14148
#   D-loop   : 1    – 576  and 16024 – 16569 (wraps around origin)
# ══════════════════════════════════════════════════════════════════════════════

het_df <- het_df %>%
  mutate(Region = case_when(
    Position >= 577   & Position <= 1601  ~ "12S_rRNA",
    Position >= 12337 & Position <= 14148 ~ "ND5",
    Position >= 16024 & Position <= 16569 ~ "D_loop",
    Position >= 1     & Position <= 576   ~ "D_loop",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Region))

cat("Heteroplasmic calls in regions of interest:", nrow(het_df), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: Per-sample detail table
# ══════════════════════════════════════════════════════════════════════════════

het_per_sample <- het_df %>%
  dplyr::select(Region, Position, REF, ALT, sample_id, 
                ancestry, macrohaplogroup, GT, AF) %>%
  arrange(Region, Position, ancestry, sample_id)

print(het_per_sample)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6: Write outputs
# ══════════════════════════════════════════════════════════════════════════════

write.csv(het_per_sample,       "het_variants_per_sample_regional.csv",  row.names = FALSE)
