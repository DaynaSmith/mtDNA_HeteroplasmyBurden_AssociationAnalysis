################################################################################################
# Negative binomial regression models to identify factors associated with heteroplasmic burden #
################################################################################################

# Load required libraries
library(MASS)  # for glm.nb
library(dplyr)

# Make sure categorical variables are factors
analysis_data$sex <- as.factor(analysis_data$sex)
analysis_data$phenotype <- as.factor(analysis_data$phenotype)
analysis_data$haplogroup <- as.factor(analysis_data$haplogroup)

# Set reference level for phenotype (controls as reference)
analysis_data$phenotype <- relevel(analysis_data$phenotype, ref = "Control")
analysis_data$haplogroup <- relevel(analysis_data$haplogroup, ref = "L0d1")

# --- Negative binomial regression for total number of heteroplasmies ---

model_total <- glm.nb(total_heteroplasmy_5_95 ~ phenotype + age + sex + haplogroup + mean_depth,
                    data = analysis_data)

# View full results
summary(model_total)

# Extract coefficients with confidence intervals
coef_summary <- summary(model_total)$coefficients

# Calculate Incidence Rate Ratios (IRR) and 95% CI
irr <- exp(coef(model_total))
irr_ci <- exp(confint(model_total))

results_total <- data.frame(
  variable = names(coef(model_total)),
  coefficient = coef(model_total),
  IRR = irr,
  CI_lower = irr_ci[,1],
  CI_upper = irr_ci[,2],
  p_value = coef_summary[,4]
)

# --- Test for age × TB_status interaction ---

model_total_interaction <- glm.nb(total_heteroplasmy_5_95 ~ phenotype * age + sex + haplogroup + mean_depth,
                                data = analysis_data)

summary(model_total_interaction)

# Extract coefficients with confidence intervals
coef_summary <- summary(model_total_interaction)$coefficients

# Calculate Incidence Rate Ratios (IRR) and 95% CI
irr <- exp(coef(model_total_interaction))
irr_ci <- exp(confint(model_total_interaction))

results_total <- data.frame(
  variable = names(coef(model_total_interaction)),
  coefficient = coef(model_total_interaction),
  IRR = irr,
  CI_lower = irr_ci[,1],
  CI_upper = irr_ci[,2],
  p_value = coef_summary[,4]
)

# Apply FDR correction (Benjamini-Hochberg)
results_total$p_adjusted_FDR <- p.adjust(results_total$p_value, method = "BH")

# --- Negative binomial regression for number of low level heteroplasmies ---

model_low <- glm.nb(low_heteroplasmy_5_10 ~ phenotype + age + sex + haplogroup + mean_depth,
                      data = analysis_data)

# View full results
summary(model_low)

# Extract coefficients with confidence intervals
coef_summary <- summary(model_low)$coefficients

# Calculate Incidence Rate Ratios (IRR) and 95% CI
irr <- exp(coef(model_low))
irr_ci <- exp(confint(model_low))

results_low <- data.frame(
  variable = names(coef(model_low)),
  coefficient = coef(model_low),
  IRR = irr,
  CI_lower = irr_ci[,1],
  CI_upper = irr_ci[,2],
  p_value = coef_summary[,4]
)

# --- Test for age × TB_status interaction ---
model_low_interaction <- glm.nb(low_heteroplasmy_5_10 ~ phenotype * age + sex + haplogroup + mean_depth,
                                  data = analysis_data)

# View full results
summary(model_low_interaction)

# Extract coefficients with confidence intervals
coef_summary <- summary(model_low_interaction)$coefficients

# Calculate Incidence Rate Ratios (IRR) and 95% CI
irr <- exp(coef(model_low_interaction))
irr_ci <- exp(confint(model_low_interaction))

results_low <- data.frame(
  variable = names(coef(model_low_interaction)),
  coefficient = coef(model_low_interaction),
  IRR = irr,
  CI_lower = irr_ci[,1],
  CI_upper = irr_ci[,2],
  p_value = coef_summary[,4]
)

# Apply FDR correction (Benjamini-Hochberg)
results_low$p_adjusted_FDR <- p.adjust(results_low$p_value, method = "BH")

# --- Negative binomial regression for number of intermediate level heteroplasmies ---

model_int <- glm.nb(intermediate_heteroplasmy_10_95 ~ phenotype + age + sex + haplogroup + mean_depth,
                    data = analysis_data)

# View full results
summary(model_int)

# Extract coefficients with confidence intervals
coef_summary <- summary(model_int)$coefficients

# Calculate Incidence Rate Ratios (IRR) and 95% CI
irr <- exp(coef(model_int))
irr_ci <- exp(confint(model_int))

results_int <- data.frame(
  variable = names(coef(model_int)),
  coefficient = coef(model_int),
  IRR = irr,
  CI_lower = irr_ci[,1],
  CI_upper = irr_ci[,2],
  p_value = coef_summary[,4]
)

# --- Test for age × TB_status interaction ---

model_int_interaction <- glm.nb(intermediate_heteroplasmy_10_95 ~ phenotype * age + sex + haplogroup + mean_depth,
                                data = analysis_data)

summary(model_int_interaction)

# View full results
summary(model_int_interaction)

# Extract coefficients with confidence intervals
coef_summary <- summary(model_int_interaction)$coefficients

# Calculate Incidence Rate Ratios (IRR) and 95% CI
irr <- exp(coef(model_int_interaction))
irr_ci <- exp(confint(model_int_interaction))

results_int <- data.frame(
  variable = names(coef(model_int_interaction)),
  coefficient = coef(model_int_interaction),
  IRR = irr,
  CI_lower = irr_ci[,1],
  CI_upper = irr_ci[,2],
  p_value = coef_summary[,4]
)

# Apply FDR correction (Benjamini-Hochberg)
results_int$p_adjusted_FDR <- p.adjust(results_int$p_value, method = "BH")

# --- Negative binomial regression for total number of homoplasmic variants ---

model_homo <- glm.nb(homoplasmy_gt95 ~ phenotype + age + sex + haplogroup + mean_depth,
                    data = analysis_data)

# View full results
summary(model_homo)

# Extract coefficients with confidence intervals
coef_summary <- summary(model_homo)$coefficients

# Calculate Incidence Rate Ratios (IRR) and 95% CI
irr <- exp(coef(model_homo))
irr_ci <- exp(confint(model_homo))

results_homo <- data.frame(
  variable = names(coef(model_homo)),
  coefficient = coef(model_homo),
  IRR = irr,
  CI_lower = irr_ci[,1],
  CI_upper = irr_ci[,2],
  p_value = coef_summary[,4]
)

# --- Test for age × TB_status interaction ---

model_homo_interaction <- glm.nb(homoplasmy_gt95 ~ phenotype * age + sex + haplogroup + mean_depth,
                                data = analysis_data)

summary(model_homo_interaction)

# Extract coefficients with confidence intervals
coef_summary <- summary(model_homo_interaction)$coefficients

# Calculate Incidence Rate Ratios (IRR) and 95% CI
irr <- exp(coef(model_homo_interaction))
irr_ci <- exp(confint(model_homo_interaction))

results_homo <- data.frame(
  variable = names(coef(model_homo_interaction)),
  coefficient = coef(model_homo_interaction),
  IRR = irr,
  CI_lower = irr_ci[,1],
  CI_upper = irr_ci[,2],
  p_value = coef_summary[,4]
)

# Apply FDR correction (Benjamini-Hochberg)
results_homo$p_adjusted_FDR <- p.adjust(results_homo$p_value, method = "BH")
