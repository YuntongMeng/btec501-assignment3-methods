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

Run sample quality checks and quantile normalization, then fit limma on the probe sets with the contrast Group 2 minus Group 1. Apply Benjamini–Hochberg correction and compare absolute signal cutoffs of 0, 0.585, and 1 at adjusted P < 0.05. These cutoffs were chosen assuming the input is on a log2 scale, where they correspond to no minimum, 1.5-fold and 2-fold changes; see [Note on the data scale](#note-on-the-data-scale). The heatmap uses the strict cutoff. For GO Biological Process ORA, map probe sets to unique Entrez IDs and use the tested, mappable genes as the background. For preranked GSEA, retain one probe per gene based on the highest mean signal and rank genes by the signed moderated t statistic.

## ExpressAnalyst differential expression check

Export the gene mapping report from the ExpressAnalyst run for the same expression matrix. Its columns are `original_id`, `entrez_id`, and `match_status` (`TRUE` or `FALSE`). Use the same mapping report as the run you want to verify.

```sh
Rscript verify_expressanalyst.R path/to/expression.tsv path/to/gene_mapping_report.csv output/expressanalyst_de.csv
```

Average the signals of matched probes with the same Entrez ID separately for each sample and retain the unmatched probe sets. Then run quantile normalization, limma with Group 2 minus Group 1, and Benjamini–Hochberg correction. The script writes the full differential expression table and prints counts for absolute signal cutoffs of 0, 0.58, and 1 at adjusted P < 0.05, the same log2-based cutoffs as the website. This check covers differential expression; the report's website plots and enrichment results come from ExpressAnalyst.

## Note on the data scale

Both workflows treat the input as log2-scaled, which is also what ExpressAnalyst assumes when it reports log2 fold changes. While reproducing the website workflow in R, I checked this assumption independently, because the values looked unusual for log2 microarray data: their range is −0.54 to 9.36, the median is 3.2, and some values are negative.

Tracing the data showed that it matches the public MAMA ColonData "Japan" matrix, derived from GEO GSE4554. That matrix is documented as log2 and was produced by MAS5 → VSN → quantile normalization. However, its numerical scale could not be reproduced exactly. Without any extra rescaling, the values are closest to a natural-log scale.

This uncertainty does not change adjusted P-values, rankings or GSEA. It changes how the fold-change cutoffs should be read, and therefore which features enter the heatmap and ORA: if the data were on a natural-log scale, a cutoff of 1 would correspond to about a 2.7-fold change rather than a 2-fold change. The full investigation, code and figures are in [`scale_audit/`](scale_audit/README.md).
