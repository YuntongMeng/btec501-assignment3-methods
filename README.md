# BTEC 501 Assignment 3: Microarray analysis methods

This repository documents two R workflows for a two-group Affymetrix Human Genome U133 Plus 2.0 microarray comparison. ExpressAnalyst was used for the assignment's main analysis. The R scripts provide an independent probe-level analysis and a separate reproduction of the website's differential expression workflow.

## Setup

Install the R packages in an R session:

```r
install.packages(c("BiocManager", "ggplot2", "pheatmap"))
BiocManager::install(c("limma", "preprocessCore", "AnnotationDbi",
                       "hgu133plus2.db", "org.Hs.eg.db", "clusterProfiler",
                       "BiocParallel", "enrichplot"))
```

Supply your own authorized input files. The expression input is a tab-separated table whose first header cell is `#NAME` and whose remaining header cells are sample IDs. The first data row starts with `#CLASS` and assigns each sample to group `1` or `2`. Each subsequent row has a unique probe-set ID followed by its numeric expression signals. The scripts operate on the supplied signal scale without an additional log transformation.

## Probe-level workflow

```sh
Rscript assignment3_checkpoint.R path/to/expression.tsv output/probe_level
```

Run sample quality checks and quantile normalization, then fit limma on the probe sets with the contrast Group 2 minus Group 1. Apply Benjamini–Hochberg correction and compare absolute signal cutoffs of 0, 0.585, and 1 at adjusted P < 0.05. The heatmap uses the strict cutoff. For GO Biological Process ORA, map probe sets to unique Entrez IDs and use the tested, mappable genes as the background. For preranked GSEA, retain one probe per gene based on the highest mean signal and rank genes by the signed moderated t statistic.

## ExpressAnalyst differential expression check

Export the gene mapping report from the ExpressAnalyst run for the same expression matrix. Its columns are `original_id`, `entrez_id`, and `match_status` (`TRUE` or `FALSE`). Use the same mapping report as the run you want to verify.

```sh
Rscript verify_expressanalyst.R path/to/expression.tsv path/to/gene_mapping_report.csv output/expressanalyst_de.csv
```

Average the signals of matched probes with the same Entrez ID separately for each sample and retain the unmatched probe sets. Then run quantile normalization, limma with Group 2 minus Group 1, and Benjamini–Hochberg correction. The script writes the full differential expression table and prints counts for absolute signal cutoffs of 0, 0.58, and 1 at adjusted P < 0.05. This check covers differential expression; the report's website plots and enrichment results come from ExpressAnalyst.

## Scale audit of the input data

The [`scale_audit/`](scale_audit/) folder traces the public source of the MAMA ColonData "Japan" matrix (GEO GSE4554) and tests which log scale it is on. See [scale_audit/README.md](scale_audit/README.md).
