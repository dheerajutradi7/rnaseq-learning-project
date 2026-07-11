#############################################
# Step 12: Import salmon quant + run DESeq2
# Project: Glioblastoma (BT935) TMZRT +/- UK5099 RNA-seq
#############################################

library(tximport)
library(DESeq2)

# ---- Paths ----
base_dir   <- "~/rnaseq_project"
quant_dir  <- file.path(base_dir, "results/salmon_quant")
tx2gene_f  <- file.path(base_dir, "reference/tx2gene_clean.tsv")
meta_f     <- file.path(base_dir, "data/sample_metadata.csv")
out_dir    <- file.path(base_dir, "results/deseq2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load metadata ----
coldata <- read.csv(meta_f, stringsAsFactors = TRUE)
rownames(coldata) <- coldata$sample
print(coldata)

# ---- Build path to each sample's quant.sf ----
files <- file.path(quant_dir, coldata$sample, "quant.sf")
names(files) <- coldata$sample
stopifnot(all(file.exists(files)))
cat("All quant.sf files found.\n")

# ---- Load tx2gene mapping ----
tx2gene <- read.table(tx2gene_f, header = FALSE, sep = "\t",
                       col.names = c("TXNAME", "GENEID"))

# ---- Import with tximport (gene-level summarization) ----
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
cat("Import complete. Gene count matrix dimensions:\n")
print(dim(txi$counts))

# ---- Build DESeq2 dataset with factorial design ----
coldata$treatment <- relevel(coldata$treatment, ref = "Naive")
coldata$drug      <- relevel(coldata$drug, ref = "DMSO")

dds <- DESeqDataSetFromTximport(txi, colData = coldata,
                                 design = ~ treatment * drug)

# ---- Pre-filter low count genes ----
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("Genes retained after filtering:", nrow(dds), "\n")

# ---- Run DESeq2 ----
dds <- DESeq(dds)

# ---- Save DESeq2 object for downstream steps ----
saveRDS(dds, file.path(out_dir, "dds.rds"))
cat("Saved DESeq2 object to", file.path(out_dir, "dds.rds"), "\n")

# ---- Results: main effect of treatment (TMZRT vs Naive) ----
res_treatment <- results(dds, name = "treatment_TMZRT_vs_Naive")
write.csv(as.data.frame(res_treatment), file.path(out_dir, "DE_treatment_TMZRT_vs_Naive.csv"))

# ---- Results: main effect of drug (UK5099 vs DMSO) ----
res_drug <- results(dds, name = "drug_UK5099_vs_DMSO")
write.csv(as.data.frame(res_drug), file.path(out_dir, "DE_drug_UK5099_vs_DMSO.csv"))

# ---- Results: interaction effect (does UK5099 change TMZRT response?) ----
res_interaction <- results(dds, name = "treatmentTMZRT.drugUK5099")
write.csv(as.data.frame(res_interaction), file.path(out_dir, "DE_interaction.csv"))

cat("\n=== Summary: treatment (TMZRT vs Naive) ===\n")
summary(res_treatment)

cat("\n=== Summary: drug (UK5099 vs DMSO) ===\n")
summary(res_drug)

cat("\n=== Summary: interaction ===\n")
summary(res_interaction)

cat("\nAll DE result tables saved in:", out_dir, "\n")
