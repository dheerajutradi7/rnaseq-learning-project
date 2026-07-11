library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

base_dir <- "~/rnaseq_project"
out_dir  <- file.path(base_dir, "results/deseq2")
enrich_dir <- file.path(base_dir, "results/enrichment")
fig_dir  <- file.path(base_dir, "figures")
dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)

run_enrichment <- function(res_csv, label) {
  cat("\n\n==============================\n")
  cat("Running enrichment for:", label, "\n")
  cat("==============================\n")

  res <- read.csv(res_csv, row.names = 1)
  res <- res[!is.na(res$padj), ]

  sig_genes <- rownames(res)[res$padj < 0.05 & abs(res$log2FoldChange) > 1]
  cat("Number of significant genes for enrichment:", length(sig_genes), "\n")

  if (length(sig_genes) < 10) {
    cat("Too few significant genes for reliable enrichment. Skipping.\n")
    return(NULL)
  }

  background <- rownames(res)

  ego <- enrichGO(gene          = sig_genes,
                   universe      = background,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)

  write.csv(as.data.frame(ego), file.path(enrich_dir, paste0("GO_BP_", label, ".csv")))

  if (nrow(as.data.frame(ego)) > 0) {
    p <- dotplot(ego, showCategory = 15) + ggtitle(paste("GO Biological Process:", label))
    ggsave(file.path(fig_dir, paste0("GO_BP_dotplot_", label, ".png")), p, width = 9, height = 7, dpi = 300)
    cat("Saved GO BP results and dotplot for", label, "\n")
  } else {
    cat("No significant GO BP terms found for", label, "\n")
  }

  gene_map <- bitr(sig_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  bg_map   <- bitr(background, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  if (nrow(gene_map) >= 10) {
    ekegg <- enrichKEGG(gene         = gene_map$ENTREZID,
                         universe     = bg_map$ENTREZID,
                         organism     = "hsa",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2)

    write.csv(as.data.frame(ekegg), file.path(enrich_dir, paste0("KEGG_", label, ".csv")))

    if (nrow(as.data.frame(ekegg)) > 0) {
      p2 <- dotplot(ekegg, showCategory = 15) + ggtitle(paste("KEGG Pathways:", label))
      ggsave(file.path(fig_dir, paste0("KEGG_dotplot_", label, ".png")), p2, width = 9, height = 7, dpi = 300)
      cat("Saved KEGG results and dotplot for", label, "\n")
    } else {
      cat("No significant KEGG pathways found for", label, "\n")
    }
  } else {
    cat("Too few mapped Entrez IDs for KEGG analysis on", label, "\n")
  }

  return(list(go = ego, kegg = if (exists("ekegg")) ekegg else NULL))
}

res_treatment_enrich <- run_enrichment(file.path(out_dir, "DE_treatment_TMZRT_vs_Naive.csv"), "treatment")
res_drug_enrich      <- run_enrichment(file.path(out_dir, "DE_drug_UK5099_vs_DMSO.csv"), "drug")
res_interaction_enrich <- run_enrichment(file.path(out_dir, "DE_interaction.csv"), "interaction")

cat("\n\nAll enrichment results saved in:", enrich_dir, "\n")
cat("All enrichment plots saved in:", fig_dir, "\n")
