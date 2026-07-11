# Bulk RNA-seq Analysis: Glioblastoma Response to Chemoradiation ± Mitochondrial Metabolism Inhibition

## Overview
This project analyzes bulk RNA-seq data from the glioblastoma cell line **BT935**, profiling gene expression across four conditions to study how blocking mitochondrial pyruvate metabolism (**UK5099**) affects the transcriptional response to standard chemoradiation therapy (**TMZRT** — Temozolomide + Radiotherapy).

## Project Nature

This is a **practice / learning project** built entirely on publicly available data (SRA BioProject PRJNA1065856). It was created to demonstrate a complete, reproducible bulk RNA-seq analysis pipeline — from raw reads to biological interpretation — as a portfolio piece. No original wet-lab data was generated; all credit for the underlying dataset belongs to the original submitters (University of Manitoba).

**Design:** 2x2 factorial, n=3 per group, 12 samples total.

| Group | Treatment | Drug |
|---|---|---|
| 1 | Naive | DMSO (vehicle) |
| 2 | Naive | UK5099 |
| 3 | TMZRT | DMSO (vehicle) |
| 4 | TMZRT | UK5099 |

**Data source:** [BioProject PRJNA1065856](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1065856) (SRA), University of Manitoba.

## Pipeline
| Step | Tool |
|---|---|
| Download | SRA prefetch / ENA direct FASTQ |
| QC | FastQC + MultiQC |
| Trimming | fastp (auto adapter detection) |
| Quantification | Salmon (selective alignment, GRCh38 cDNA, Ensembl release 110) |
| Gene-level import | tximport |
| Differential expression | DESeq2, factorial design `~ treatment * drug` |
| Enrichment | clusterProfiler (GO Biological Process, KEGG) |
| Visualization | ggplot2, pheatmap |

## Key Results

- 12/12 samples passed QC with high quality (Phred >28 throughout), minimal adapter content.
- 90-92% mapping rate for all samples (Salmon selective alignment).
- 16,141 genes retained after low-count filtering (from 38,366 total).

| Comparison | Significant genes (padj<0.1) | % of tested genes |
|---|---|---|
| TMZRT vs Naive | 9,126 | 57% |
| UK5099 vs DMSO | 5,040 | 31% |
| Interaction (TMZRT x UK5099) | 1,564 | 9.7% |

- PCA shows clean separation of all four groups, with treatment explaining 70% of variance (PC1) and drug condition explaining 9% (PC2).
- GO/KEGG enrichment for the TMZRT effect is dominated by immune/inflammatory and cytokine-signaling pathways.
- The UK5099 effect enriches for angiogenesis, epithelial differentiation, and extracellular matrix pathways.
- The interaction analysis shows UK5099 measurably changes the transcriptional response to TMZRT.

## Data Availability
Raw and intermediate files (FASTQ, salmon index, quant.sf) are not included in this repository due to size. They can be regenerated using the accession list in `data/sample_metadata.csv` and the pipeline steps above, or downloaded directly from SRA BioProject PRJNA1065856.

## Environment
```bash
conda create -n rnaseq -c bioconda -c conda-forge sra-tools fastqc multiqc fastp salmon
conda create -n r_deseq -c bioconda -c conda-forge r-base=4.3 bioconductor-deseq2 bioconductor-tximport r-ggplot2 r-pheatmap bioconductor-clusterprofiler bioconductor-org.hs.eg.db
```

## Walkthrough: Reproducing This Pipeline With Your Own Samples

The steps below are general enough to apply to any paired-end bulk RNA-seq dataset, not just this one.

### 1. Set up environments
```bash
conda create -n rnaseq -c bioconda -c conda-forge sra-tools fastqc multiqc fastp salmon
conda create -n r_deseq -c bioconda -c conda-forge r-base=4.3 bioconductor-deseq2 bioconductor-tximport r-ggplot2 r-pheatmap bioconductor-clusterprofiler bioconductor-org.hs.eg.db
```

### 2. Get your FASTQ files
- From SRA: `prefetch --option-file your_srr_list.txt --output-directory ./raw_sra`, then `fasterq-dump` to convert to FASTQ
- Or from ENA directly (often faster): `https://ftp.sra.ebi.ac.uk/vol1/fastq/<prefix>/<SRR>/<SRR>_1.fastq.gz`
- Or your own sequencing facility's FASTQ output — just make sure files are named consistently as `SAMPLE_1.fastq.gz` / `SAMPLE_2.fastq.gz` (paired-end)

### 3. Quality check
```bash
fastqc data/raw_fastq/*.fastq.gz -o results/fastqc_raw -t 4
multiqc results/fastqc_raw -o results/multiqc_raw
```
Review the MultiQC report — check adapter content, quality scores, and GC distribution before deciding on trimming stringency.

### 4. Trim adapters/low quality
```bash
fastp -i R1.fastq.gz -I R2.fastq.gz -o R1.trimmed.fastq.gz -O R2.trimmed.fastq.gz \
  --detect_adapter_for_pe --thread 4 --json sample.json --html sample.html
```
Loop this over all sample pairs. Re-run FastQC/MultiQC on trimmed output to confirm cleanup.

### 5. Build a reference index (once, matching your organism)
```bash
wget <Ensembl cDNA FASTA for your organism>
salmon index -t reference.cdna.fa.gz -i salmon_index -k 31 -p 4
```

### 6. Quantify each sample
```bash
salmon quant -i salmon_index -l A -1 R1.trimmed.fastq.gz -2 R2.trimmed.fastq.gz -p 4 -o results/salmon_quant/SAMPLE
```

### 7. Build a tx2gene map (matches your reference)
```bash
zcat reference.cdna.fa.gz | grep "^>" | awk '{print $1"\t"$4}' | sed 's/>//; s/gene://; s/gene_symbol://; s/\.[0-9]*//g' > tx2gene.tsv
```

### 8. Prepare sample metadata
Create `sample_metadata.csv` with one row per sample, columns for `sample` and your experimental factors (e.g., `treatment`, `drug`, `timepoint`) — this defines your DESeq2 design formula.

### 9. Run differential expression (R)
Use `scripts/01_import_deseq2.R` as a template — update `design = ~ your_factors`, sample metadata path, and comparison names to match your experiment.

### 10. Visualize and interpret
Use `scripts/02_visualization.R` (PCA, volcano, heatmap) and `scripts/03_enrichment.R` (GO/KEGG) as templates — both just need the DE result CSV paths updated to match your comparisons.

### Notes for adapting to your data
- Adjust the DESeq2 `design` formula to match your actual experimental factors (this project used a 2x2 factorial; a simple two-group study would use `~ condition`)
- For non-human organisms, swap `org.Hs.eg.db` for the appropriate Bioconductor annotation package (e.g., `org.Mm.eg.db` for mouse) and use `organism = "mmu"` in `enrichKEGG`
- Sample size affects statistical power — this project used n=3/group; more replicates generally give more reliable results
