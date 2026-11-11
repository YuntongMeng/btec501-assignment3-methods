# 02  Provenance: does the MAMA Japan matrix come from GSE4554, and how?
#   - group labels in MAMA vs. GEO sample descriptions
#   - per-sample rank agreement with the original MAS5 signal (Spearman)
#   - equal MAS5 values within a sample -> equal MAMA values? (a per-sample function)
#   - negative MAMA values vs. MAS5 detection calls (Absent / Present)
source("scale_audit/R/common.R")
suppressPackageStartupMessages(library(GEOquery))
x <- load_inputs()
J <- x$japan; smp <- colnames(J); p <- rownames(J)
M <- x$mas5[p, smp]

geo_status <- sub(",.*", "", x$pdata[smp, "description"])
geo_site   <- trimws(sub(".*,", "", x$pdata[smp, "description"]))

# Detection calls are only in the per-sample GSM tables
calls <- sapply(smp, function(s) {
  t <- Table(getGEO(s, destdir = DATA))
  setNames(as.character(t$ABS_CALL), t$ID_REF)[p]
})

per_sample <- data.frame(
  GSM = smp,
  MAMA_status = as.character(x$japan_cl[smp, "MSIstatus"]),
  GEO_status = geo_status, GEO_site = geo_site,
  spearman_vs_MAS5 = sapply(smp, function(s) cor(J[, s], M[, s], method = "spearman")),
  max_spread_equal_MAS5 = sapply(smp, function(s)          # 0 => a function of MAS5
    max(tapply(J[, s], M[, s], function(v) diff(range(v))))),
  n_negative = colSums(J < 0),
  row.names = NULL)
write.csv(per_sample, file.path(RESULTS, "provenance_per_sample.csv"), row.names = FALSE)

neg <- J < 0
summary_tab <- data.frame(
  check = c("Samples in MAMA Japan", "Labels agreeing with GEO (MSI/MSS)",
            "Min per-sample Spearman vs MAS5", "Max spread among equal MAS5 values",
            "Negative values", "Negative values with MAS5 call = Absent",
            "MAS5 signal behind negative values (median)",
            "MAS5 signal behind negative values (max)"),
  value = c(length(smp), sum(per_sample$MAMA_status == per_sample$GEO_status),
            round(min(per_sample$spearman_vs_MAS5), 6),
            max(per_sample$max_spread_equal_MAS5),
            sum(neg), sum(calls[neg] == "A"), median(M[neg]), max(M[neg])))
write.csv(summary_tab, file.path(RESULTS, "provenance_summary.csv"), row.names = FALSE)
print(summary_tab, row.names = FALSE)
print(table(MAMA = per_sample$MAMA_status, GEO = per_sample$GEO_status))
