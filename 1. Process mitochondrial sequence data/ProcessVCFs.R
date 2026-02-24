library(vcfR)
library(tidyverse)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: Load data
# ══════════════════════════════════════════════════════════════════════════════

vcf            <- read.vcfR("LTBI-ResisTB_merged_MaxMissing0.95_PASS_FinalSampleList_OnlyAfricanSamples.recode.vcf")
haplo_defining <- read.table("per_sample_haplodefining_positions.tsv",
                             header = TRUE, sep = "\t")

vcf_positions  <- as.integer(vcf@fix[, "POS"])
gt_matrix      <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
vcf_samples    <- colnames(gt_matrix)

# Sanity check sample name matching
missing_from_haplo <- setdiff(vcf_samples, unique(haplo_defining$SampleID))
missing_from_vcf   <- setdiff(unique(haplo_defining$SampleID), vcf_samples)
if (length(missing_from_haplo) > 0) warning("In VCF but not haplo file: ", paste(missing_from_haplo, collapse = ", "))
if (length(missing_from_vcf)   > 0) warning("In haplo file but not VCF: ", paste(missing_from_vcf,   collapse = ", "))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Build per-sample mask of haplogroup-defining homoplasmic calls
#   Homoplasmic alt = GT "1" or "2"  →  mask these if haplogroup-defining
#   Heteroplasmic   = GT "0/1"/"0/2" →  always preserved
# ══════════════════════════════════════════════════════════════════════════════

haplo_mask <- matrix(FALSE, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix),
                     dimnames = list(NULL, vcf_samples))

for (samp in vcf_samples) {
  haplo_pos     <- haplo_defining$Position[haplo_defining$SampleID == samp]
  rows_to_check <- which(vcf_positions %in% haplo_pos)
  
  if (length(rows_to_check) > 0) {
    gt_vals            <- gt_matrix[rows_to_check, samp]
    is_homoplasmic_alt <- !is.na(gt_vals) & gt_vals %in% c("1", "2")
    haplo_mask[rows_to_check[is_homoplasmic_alt], samp] <- TRUE
  }
}

cat("Haplogroup-defining homoplasmic calls to mask:", sum(haplo_mask), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: Apply mask to raw @gt slot, zeroing GT and AF for masked entries
#         DP is preserved. FORMAT field ordering is handled per-row since
#         it is inconsistent (GT:DP:AF vs GT:AF:DP) across rows.
# ══════════════════════════════════════════════════════════════════════════════

raw_gt <- vcf@gt

for (samp in vcf_samples) {
  col_idx      <- which(colnames(raw_gt) == samp)
  rows_to_mask <- which(haplo_mask[, samp])
  
  if (length(rows_to_mask) > 0) {
    for (row_idx in rows_to_mask) {
      entry <- raw_gt[row_idx, col_idx]
      if (is.na(entry)) next
      
      format_fields <- strsplit(raw_gt[row_idx, 1], ":")[[1]]
      entry_fields  <- strsplit(entry, ":")[[1]]
      
      gt_idx <- which(format_fields == "GT")
      af_idx <- which(format_fields == "AF")
      
      if (length(gt_idx) > 0) entry_fields[gt_idx] <- "0"   # set to ref
      if (length(af_idx) > 0) entry_fields[af_idx] <- "0"   # set AF to 0
      
      raw_gt[row_idx, col_idx] <- paste(entry_fields, collapse = ":")
    }
  }
}

vcf_filtered     <- vcf
vcf_filtered@gt  <- raw_gt

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: Verify heteroplasmic calls are preserved
# ══════════════════════════════════════════════════════════════════════════════

gt_filtered <- extract.gt(vcf_filtered, element = "GT", as.numeric = FALSE)

het_before <- sum(gt_matrix  %in% c("0/1", "0/2"), na.rm = TRUE)
het_after  <- sum(gt_filtered %in% c("0/1", "0/2"), na.rm = TRUE)
cat("Heteroplasmic calls before:", het_before, "\n")
cat("Heteroplasmic calls after: ", het_after,  "\n")
if (het_before != het_after) warning("Heteroplasmic call count changed — investigate before proceeding!")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: Extract AF and classify variants per sample
#
#   Using AF for classification since it is reliably populated for all
#   non-reference calls and allows heteroplasmy subcategories.
#
#   Categories:
#     reference              AF == 0   (not a variant)
#     very_low_heteroplasmy  0 < AF < 0.05
#     low_heteroplasmy       0.05 <= AF < 0.10
#     intermediate_hetero    0.10 <= AF < 0.95
#     homoplasmy             AF >= 0.95
# ══════════════════════════════════════════════════════════════════════════════

af_matrix  <- extract.gt(vcf_filtered, element = "AF", as.numeric = FALSE)

# Convert to numeric, treating "." and NA as 0 (reference)
clean_af <- function(x) {
  if (is.na(x) || x == "." || x == "") return(0)
  as.numeric(x)
}
vaf_matrix <- apply(af_matrix, c(1, 2), clean_af)

# Categorise
categorize_vaf <- function(vaf) {
  case_when(
    vaf == 0                  ~ "reference",
    vaf > 0   & vaf < 0.05   ~ "very_low_heteroplasmy",
    vaf >= 0.05 & vaf < 0.10 ~ "low_heteroplasmy",
    vaf >= 0.10 & vaf < 0.95 ~ "intermediate_heteroplasmy",
    vaf >= 0.95               ~ "homoplasmy",
    TRUE                      ~ "other"
  )
}

vaf_category_matrix <- apply(vaf_matrix, c(1, 2), categorize_vaf)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6: Count variants per category per sample
# ══════════════════════════════════════════════════════════════════════════════

variant_counts <- data.frame(
  sample_id                       = vcf_samples,
  homoplasmy_gt95                 = colSums(vaf_category_matrix == "homoplasmy",                 na.rm = TRUE),
  intermediate_heteroplasmy_10_95 = colSums(vaf_category_matrix == "intermediate_heteroplasmy",  na.rm = TRUE),
  low_heteroplasmy_5_10           = colSums(vaf_category_matrix == "low_heteroplasmy",           na.rm = TRUE),
  total_heteroplasmy_5_95         = colSums(vaf_category_matrix == "low_heteroplasmy",           na.rm = TRUE) +
    colSums(vaf_category_matrix == "intermediate_heteroplasmy",  na.rm = TRUE),
  very_low_heteroplasmy_lt5       = colSums(vaf_category_matrix == "very_low_heteroplasmy",      na.rm = TRUE)
)

cat("\nVariant count summary after haplogroup-defining site removal:\n")
cat("Mean homoplasmic variants per sample:           ", round(mean(variant_counts$homoplasmy_gt95),                 2), "\n")
cat("Mean heteroplasmic variants per sample (5-95%): ", round(mean(variant_counts$total_heteroplasmy_5_95),         2), "\n")
cat("Mean very low VAF variants per sample (<5%):    ", round(mean(variant_counts$very_low_heteroplasmy_lt5),       2), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7: Write outputs
# ══════════════════════════════════════════════════════════════════════════════

write.csv(variant_counts, "VariantCountsByVAF.csv", row.names = FALSE)
write.vcf(vcf_filtered,   "vcf_haplo_filtered.vcf.gz")

cat("\nDone.\n")
cat("  Variant counts: VariantCountsByVAF.csv\n")
cat("  Filtered VCF:   vcf_haplo_filtered.vcf.gz\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 8: Merge with metadata 
# ══════════════════════════════════════════════════════════════════════════════

sample_data <- read.csv("Heteroplasmy_analysis/Merged_cohorts/sample_data.csv", header=TRUE)
analysis_data <- merge(sample_data, variant_counts, by="sample_id")

head(analysis_data)
