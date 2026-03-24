# Script: 06_analysis.R
# Purpose: Diversity analysis and visualization of gut microbiome data
# Input: combined_bracken.txt
# Output: figures/ directory with all plots
# Author: Vian Lelo | Date: March 23rd 2026
# NOTE: Each figure will display in RStudio AND be saved to figures/

# 0. PACKAGE INSTALLATION (RUN ONCE, NOT EVERY TIME)

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("ANCOMBC", "microbiome", "phyloseq"))
# Optional: if you specifically need this CVXR version
# remotes::install_version("CVXR", version = "1.0-11")


# 1. LOAD LIBRARIES

library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(ANCOMBC)
library(microbiome)


# 2. SETUP OUTPUT DIRECTORIES

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)


# 3. LOAD DATA

data <- read.table("combined_bracken.txt",
                   header = TRUE, sep = "\t", row.names = 1)

count_cols <- grep("_num$", colnames(data))
counts <- data[, count_cols]

colnames(counts) <- gsub(".bracken_num", "", colnames(counts))

metadata <- data.frame(
  sample = c("SRR8146951", "SRR8146952", "SRR8146944",
             "SRR8146936", "SRR8146935", "SRR8146938"),
  diet   = c("vegan", "vegan", "vegan",
             "omnivore", "omnivore", "omnivore")
)
rownames(metadata) <- metadata$sample


# 4. CREATE PHYLOSEQ OBJECT

OTU  <- otu_table(as.matrix(counts), taxa_are_rows = TRUE)
SAMP <- sample_data(metadata)
ps   <- phyloseq(OTU, SAMP)

cat("Phyloseq object created!\n")
print(ps)


# 5. TAXONOMIC ABUNDANCE BAR PLOT

top20 <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:20])
ps_top20 <- prune_taxa(top20, ps)

ps_top20_rel <- transform_sample_counts(ps_top20, function(x) x / sum(x))

fig1 <- plot_bar(ps_top20_rel, x = "sample", fill = "OTU") +
  facet_wrap(~diet, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.text = element_text(size = 7)) +
  labs(title = "Top 20 Taxa Relative Abundance by Diet Group",
       x = "Sample", y = "Relative Abundance",
       fill = "Species")

print(fig1)
ggsave("figures/01_taxonomic_abundance.pdf", fig1, width = 12, height = 7)
ggsave("figures/01_taxonomic_abundance.png", fig1, width = 12, height = 7, dpi = 300)
cat("Figure 1 saved: 01_taxonomic_abundance.pdf + .png\n")


# 6. ALPHA DIVERSITY

fig2 <- plot_richness(ps,
                      x = "diet",
                      measures = c("Observed", "Shannon", "Simpson")) +
  geom_boxplot(aes(fill = diet), alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2) +
  scale_fill_manual(values = c("vegan" = "#2C5F2D", "omnivore" = "#B85042")) +
  theme_minimal() +
  labs(title = "Alpha Diversity by Diet Group",
       x = "Diet", y = "Diversity Measure",
       fill = "Diet")

print(fig2)
ggsave("figures/02_alpha_diversity.pdf", fig2, width = 8, height = 6)
ggsave("figures/02_alpha_diversity.png", fig2, width = 8, height = 6, dpi = 300)
cat("Figure 2 saved: 02_alpha_diversity.pdf + .png\n")


# 7. BETA DIVERSITY (PCoA)

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

ord <- ordinate(ps_rel, method = "PCoA", distance = "bray")

fig3 <- plot_ordination(ps_rel, ord, color = "diet") +
  geom_point(size = 5) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("vegan" = "#2C5F2D", "omnivore" = "#B85042")) +
  theme_minimal() +
  labs(title = "Beta Diversity: Bray-Curtis PCoA",
       subtitle = "Vegan vs Omnivore Gut Microbiome",
       color = "Diet")

print(fig3)
ggsave("figures/03_beta_diversity_PCoA.pdf", fig3, width = 8, height = 6)
ggsave("figures/03_beta_diversity_PCoA.png", fig3, width = 8, height = 6, dpi = 300)
cat("Figure 3 saved: 03_beta_diversity_PCoA.pdf + .png\n")

metadata <- data.frame(sample_data(ps))
dist_matrix <- phyloseq::distance(ps_rel, method = "bray")
diet_vector <- metadata$diet
permanova <- adonis2(dist_matrix ~ diet_vector, permutations = 999)
print(permanova)
                                  
# 8. DIFFERENTIAL ABUNDANCE (ANCOM-BC2)

ancom_out <- ancombc2(
  data = ps,
  fix_formula = "diet",
  group = "diet",
  p_adj_method = "BH",
  prv_cut = 0.10,
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  verbose = TRUE
)

res_ancom <- ancom_out$res

res_sig <- res_ancom[res_ancom$diff_dietvegan == TRUE, ]
cat("Significant taxa (FDR < 0.05):", nrow(res_sig), "\n")

if (nrow(res_sig) == 0) {
  cat("No significant taxa at FDR < 0.05 — showing top 20 by p-value (exploratory)\n")
  res_sig <- res_ancom %>%
    arrange(p_dietvegan) %>%
    head(20) %>%
    mutate(direction = ifelse(lfc_dietvegan > 0,
                              "Higher in Vegan",
                              "Higher in Omnivore"))
  subtitle_text <- "ANCOM-BC2 | Top 20 by p-value | exploratory n=3"
} else {
  res_sig <- res_sig %>%
    mutate(direction = ifelse(lfc_dietvegan > 0,
                              "Higher in Vegan",
                              "Higher in Omnivore"))
  subtitle_text <- paste("ANCOM-BC2 | FDR < 0.05 |", nrow(res_sig), "taxa")
}

fig4 <- ggplot(res_sig,
               aes(x = reorder(taxon, lfc_dietvegan),
                   y = lfc_dietvegan,
                   fill = direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Higher in Vegan"    = "#2C5F2D",
                               "Higher in Omnivore" = "#B85042")) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom") +
  labs(title    = "Differential Abundance: Vegan vs Omnivore",
       subtitle = subtitle_text,
       x        = "Species",
       y        = "Log Fold Change (Vegan vs Omnivore)",
       fill     = "Direction")

print(fig4)
ggsave("figures/04_differential_abundance.pdf", fig4, width = 10, height = 8)
ggsave("figures/04_differential_abundance.png", fig4, width = 10, height = 8, dpi = 300)
cat("Figure 4 saved!\n")

write.table(res_ancom,
            "results/ancombc2_results.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Results saved: results/ancombc2_results.txt\n")
