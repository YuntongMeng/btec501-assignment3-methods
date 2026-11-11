# Scale audit: MAMA ColonData "Japan" (GSE4554)

This folder checks the provenance and numerical scale of the `Japan` component of the `ColonData` example in the R package MAMA. That matrix is documented as log2-transformed and derives from GEO GSE4554 (MAS5 signal).

**Findings in brief:** the source and the MSI/MSS labels are confirmed. The values are a per-sample monotone transform of MAS5, consistent with the documented VSN + quantile-normalization pipeline. Their numerical scale, however, cannot be reproduced exactly: they are closest to a natural-log scale, while the documented glog2 would need an unexplained rescaling. Significance and rankings are unaffected; fold-change cutoffs are not.

- Full report: [REPORT.md](REPORT.md)
- Run everything from the repository root: `Rscript scale_audit/R/run_all.R`
- Inputs are downloaded at run time from the CRAN archive and GEO into `scale_audit/data/`, which is not tracked.
- Outputs: `results/` (CSV tables, `sessionInfo.txt`) and `figures/` (PNG).

Required packages: `GEOquery`, `vsn`, `preprocessCore`, `limma` (Bioconductor) and `ggplot2`.

## Use of AI tools

I identified the scale issue myself during independent verification; AI tools were used to assist in tracing the GEO metadata, locating the MAMA package source, verifying the original author's thesis, and drafting code and text, and all results were run and checked by me.
