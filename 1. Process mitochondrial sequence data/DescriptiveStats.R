library(patchwork)
library(ggplot2)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Heteroplasmy by cohort
# ══════════════════════════════════════════════════════════════════════════════

wilcoxon_count <- wilcox.test(total_heteroplasmy_5_95 ~ cohort, 
                              data = analysis_data)
print(wilcoxon_count)

p1 <- ggplot(analysis_data, aes(x = cohort, y = total_heteroplasmy_5_95, 
                                     fill = cohort)) +
  geom_boxplot(width = 0.3, alpha = 1) +
  labs(subtitle = paste0("Wilcoxon rank-sum p = 0.003188"),
       x = "Cohort",
       y = "Number of Heteroplasmies"
  ) +
  scale_fill_manual(values = c("skyblue", "sienna1")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
p1


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Heteroplasmy by sex
# ══════════════════════════════════════════════════════════════════════════════

wilcoxon_count <- wilcox.test(total_heteroplasmy_5_95 ~ sex, 
                              data = analysis_data)
print(wilcoxon_count)

p2 <- ggplot(analysis_data, aes(x = sex, y = total_heteroplasmy_5_95, 
                                     fill = sex)) +
  geom_boxplot(width = 0.3, alpha = 1) +
  labs(subtitle = paste0("Wilcoxon rank-sum p = 0.3121"),
       x = "Cohort",
       y = "Number of Heteroplasmies"
  ) +
  scale_fill_manual(values = c("lightpink", "lightblue")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
p2

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: Heteroplasmy by age
# ══════════════════════════════════════════════════════════════════════════════

# Age correlation
age_cor <- cor.test(analysis_data$age,
                    analysis_data$total_heteroplasmy_5_95,
                    method = "spearman")

cat("Age correlation:\n")
cat("rho =", round(age_cor$estimate, 3), "\n")
cat("p-value =", format.pval(age_cor$p.value, digits = 3), "\n\n")

p3 <- ggplot(analysis_data, aes(x = age, y = total_heteroplasmy_5_95)) +
  geom_point(alpha = 0.5, color = "darkgreen") +
  geom_smooth(method = "lm", se = TRUE, color = "grey") +
  labs(
    subtitle = paste0("Spearman ρ = ", round(age_cor$estimate, 3),
                      ", p = ", format.pval(age_cor$p.value, digits = 3)),
    x = "Age (years)",
    y = "Number of Heteroplasmies"
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"))
p3

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: Heteroplasmy by sequencing depth
# ══════════════════════════════════════════════════════════════════════════════

# Read depth correlation
depth_cor <- cor.test(analysis_data$mean_read_depth,
                      analysis_data$total_heteroplasmy_5_95,
                      method = "spearman")

cat("Depth correlation:\n")
cat("rho =", round(depth_cor$estimate, 3), "\n")
cat("p-value =", format.pval(depth_cor$p.value, digits = 3), "\n\n")

p4 <- ggplot(analysis_data, aes(x = mean_read_depth, y = total_heteroplasmy_5_95)) +
  geom_point(alpha = 0.5, color = "midnightblue") +
  geom_smooth(method = "lm", se = TRUE, color = "grey") +
  labs(
    subtitle = paste0("Spearman ρ = ", round(depth_cor$estimate, 3),
                      ", p = ", format.pval(depth_cor$p.value, digits = 3)),
    x = "Mean read depth",
    y = "Number of Heteroplasmies"
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"))
p4

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5: Heteroplasmy by ancestry
# ══════════════════════════════════════════════════════════════════════════════

# Statistical test
wilcoxon_ancestry <- wilcox.test(total_heteroplasmy_5_95 ~ ancestry, 
                                 data = analysis_data)
cat("\nWilcoxon rank-sum test: p =", format.pval(wilcoxon_ancestry$p.value, digits = 3), "\n\n")

p5 <- ggplot(analysis_data, aes(x = ancestry, y = total_heteroplasmy_5_95, 
                                     fill = ancestry)) +
  geom_boxplot(width = 0.3, alpha = 1) +
  labs(
    subtitle = paste0("Wilcoxon rank-sum p = ", 
                      format.pval(wilcoxon_ancestry$p.value, digits = 3)),
    x = "Ancestry",
    y = "Number of Heteroplasmies"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
p5

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 6: Heteroplasmy by haplogroup
# ══════════════════════════════════════════════════════════════════════════════

kw_haplogroup <- kruskal.test(total_heteroplasmy_5_95 ~ macrohaplogroup,
                              data = analysis_data)
cat("\nKruskal-Wallis test: p =", format.pval(kw_haplogroup$p.value, digits = 3), "\n\n")

haplo_colors <- c(
  "L0a1"="#1B9E77",
  "L0a2"="#66C2A5",
  "L0d1"="#D95F02",
  "L0d2"="#FC8D62",
  "L0d3"="#B23A00",
  "L0f1"="#7570B3",
  "L1"="#E7298A",
  "L2"="#66A61E",
  "L3"="#E6AB02",
  "L4"="#A6761D",
  "L5"="#666666"
)


p6 <- analysis_data %>%
  ggplot(aes(x = macrohaplogroup,
             y = total_heteroplasmy_5_95, fill = macrohaplogroup)) +
  geom_boxplot(width = 0.3, alpha = 0.8) +
  labs(
    subtitle = paste0("Kruskal-Wallis p = ", 
                      format.pval(kw_haplogroup$p.value, digits = 3)),
    x = "Haplogroup",
    y = "Number of Heteroplasmies"
  ) +
  scale_fill_manual(values = haplo_colors) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) 
p6


(p1 | p2) / (p3 | p4) / (p5 | p6)
ggsave("Heteroplasmy_stats.png",
       (p1 | p2) / (p3 | p4) / (p5 | p6), width = 12, height = 12, dpi = 300)
