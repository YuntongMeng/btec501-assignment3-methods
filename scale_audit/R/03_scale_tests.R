# 03  Test three scale hypotheses against the MAMA Japan values.
#   H_ln    : natural log of MAS5, then quantile normalisation
#   H_log2  : log2 of MAS5, then quantile normalisation
#   H_glog2 : VSN (glog2, vsn::vsn2) of MAS5, then quantile normalisation
#             (the pipeline described in the MAMA author's thesis, section 5.2)
# Each candidate is computed on all 84 GSE4554 arrays, as in the original study.
#
# R^2 of a free linear fit measures SHAPE only (y = a*x + b for any a, b).
# Agreement in VALUE needs slope ~ 1, intercept ~ 0 and small absolute differences.
source("scale_audit/R/common.R")
suppressPackageStartupMessages({ library(vsn); library(limma) })
x <- load_inputs()
J <- x$japan; smp <- colnames(J); p <- rownames(J)
E <- x$mas5
floor_pos <- min(E[E > 0])                 # MAS5 contains a few zeros; avoid log(0)

cand <- list(
  ln    = qn(log(pmax(E, floor_pos))),
  log2  = qn(log2(pmax(E, floor_pos))),
  glog2 = qn(exprs(vsn2(E, verbose = FALSE))))
saveRDS(lapply(cand, function(m) m[p, smp]), file.path(DATA, "candidates.rds"))

# Per-sample fits -------------------------------------------------------------
fits <- do.call(rbind, lapply(names(cand), function(h) {
  P <- cand[[h]][p, smp]
  do.call(rbind, lapply(smp, function(s) {
    f <- lm(J[, s] ~ P[, s])
    data.frame(hypothesis = h, GSM = s, slope = coef(f)[[2]], intercept = coef(f)[[1]],
               R2 = summary(f)$r.squared, median_abs_diff = median(abs(J[, s] - P[, s])))
  }))
}))
write.csv(fits, file.path(RESULTS, "scale_fits_per_sample.csv"), row.names = FALSE)

# Level and spread checks (no fitting) ----------------------------------------
lvl <- data.frame(
  hypothesis = c("MAMA Japan (observed)", names(cand)),
  median = c(median(J), sapply(cand, function(m) median(m[p, smp]))),
  SD     = c(sd(J),     sapply(cand, function(m) sd(m[p, smp]))))

# Group-difference check: MSS - MSI estimated on each scale --------------------
grp <- factor(x$japan_cl[smp, "MSIstatus"], levels = c("MSI", "MSS"))
d <- model.matrix(~grp)
lfc <- function(m) topTable(eBayes(lmFit(m, d)), coef = 2, number = Inf, sort.by = "none")$logFC
y <- lfc(J)
fc_slope <- sapply(names(cand), function(h) coef(lm(y ~ lfc(cand[[h]][p, smp])))[[2]])

summ <- do.call(rbind, lapply(names(cand), function(h) {
  f <- fits[fits$hypothesis == h, ]
  data.frame(hypothesis = h,
             R2_median = median(f$R2), R2_min = min(f$R2),
             slope_median = median(f$slope), slope_range = paste(sprintf("%.3f", range(f$slope)), collapse = "-"),
             intercept_median = median(f$intercept),
             median_abs_diff = median(f$median_abs_diff),
             group_difference_slope = fc_slope[[h]])
}))
write.csv(summ, file.path(RESULTS, "scale_summary.csv"), row.names = FALSE)
write.csv(lvl,  file.path(RESULTS, "scale_level_spread.csv"), row.names = FALSE)
print(summ, row.names = FALSE, digits = 4)
print(lvl, row.names = FALSE, digits = 4)
cat(sprintf("Median back-transform: exp(%.3f) = %.1f vs MAS5 median %.1f (2^ = %.1f)\n",
            median(J), exp(median(J)), median(E[p, smp]), 2^median(J)))
