#############################################
# Step 13: Visualization
# PCA, Volcano plots, Heatmap
#############################################

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)

base_dir <- "~/rnaseq_project"
out_dir  <- file.path(base_dir, "results/deseq2")
fig_dir  <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load saved DESeq2 object ----
dds <- readRDS(file.path(out_dir, "dds.rds"))

# ---- Variance-stabilizing transform (for PCA/heatmap) ----
vsd <- vst(dds, blind = TRUE)

#############################################
# 1. PCA plot
#############################################
pca_data <- plotPCA(vsd, intgroup = c("treatment", "drug"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = treatment, shape = drug)) +
  geom_point(size = 4, alpha = 0.9) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of samples (VST-transformed counts)") +
  theme_bw(base_size = 13)

ggsave(file.path(fig_dir, "PCA_plot.png"), p_pca, width = 7, height = 5.5, dpi = 300)
cat("Saved PCA_plot.png\n")

#############################################
# 2. Volcano plots (one per comparison)
#############################################
make_volcano <- function(res_csv, title, out_png) {
  res <- read.csv(res_csv, row.names = 1)
  res <- res[!is.na(res$padj), ]
  res$sig <- "Not significant"
  res$sig[res$padj < 0.05 & res$log2FoldChange > 1]  <- "Up"
  res$sig[res$padj < 0.05 & res$log2FoldChange < -1] <- "Down"

  p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "Not significant" = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    theme_bw(base_size = 13) +
    ggtitle(title) +
    xlab("log2 Fold Change") + ylab("-log10 adjusted p-value")

  ggsave(out_png, p, width = 7, height = 6, dpi = 300)
  cat("Saved", out_png, "\n")
  return(res)
}

res_treat <- make_volcano(file.path(out_dir, "DE_treatment_TMZRT_vs_Naive.csv"),
                           "Volcano: TMZRT vs Naive",
                           file.path(fig_dir, "volcano_treatment.png"))

res_drug <- make_volcano(file.path(out_dir, "DE_drug_UK5099_vs_DMSO.csv"),
                          "Volcano: UK5099 vs DMSO",
                          file.path(fig_dir, "volcano_drug.png"))

res_inter <- make_volcano(file.path(out_dir, "DE_interaction.csv"),
                           "Volcano: Interaction (TMZRT x UK5099)",
                           file.path(fig_dir, "volcano_interaction.png"))

#############################################
# 3. Heatmap of top 30 DE genes (treatment comparison)
#############################################
res_treat_ordered <- res_treat[order(res_treat$padj), ]
top_genes <- rownames(res_treat_ordered)[1:30]

mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)  # center per gene

annotation_col <- as.data.frame(colData(vsd)[, c("treatment", "drug")])

png(file.path(fig_dir, "heatmap_top30_treatment.png"), width = 1600, height = 2000, res = 200)
pheatmap(mat,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 7,
         main = "Top 30 DE genes: TMZRT vs Naive")
dev.off()
cat("Saved heatmap_top30_treatment.png\n")

cat("\nAll figures saved in:", fig_dir, "\n")
