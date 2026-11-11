# 04  How much does the scale question matter for a differential expression analysis?
#   limma, MSS vs MSI, on the MAMA Japan values (41 samples).
#   - adjusted P-values and rankings are invariant to multiplying the data by a constant
#   - only fold-change cutoffs move: a log2 cutoff c corresponds to c * ln(2) on a ln scale
#   Also compares running an extra quantile normalisation (as web tools do by default).
source("scale_audit/R/common.R")
suppressPackageStartupMessages(library(limma))
x <- load_inputs()
J <- x$japan; smp <- colnames(J)
grp <- factor(x$japan_cl[smp, "MSIstatus"], levels = c("MSI", "MSS"))
d <- model.matrix(~grp)
tt <- function(m) topTable(eBayes(lmFit(m, d)), coef = 2, number = Inf, sort.by = "none")

res <- list(as_is = tt(J), extra_QN = tt(qn(J)))
r2  <- tt(J * (1 / log(2)))               # same data rescaled ln -> log2
inv <- c(max_diff_adjP_after_rescale = max(abs(res$as_is$adj.P.Val - r2$adj.P.Val)),
         cor_logFC_extraQN = cor(res$as_is$logFC, res$extra_QN$logFC))

cuts <- c(0, 0.585, 1)
sens <- do.call(rbind, lapply(names(res), function(k) {
  r <- res[[k]]
  do.call(rbind, lapply(cuts, function(cc) data.frame(
    normalisation = k, log2FC_cutoff = cc,
    n_if_scale_is_log2 = sum(r$adj.P.Val < 0.05 & abs(r$logFC) >= cc),
    n_if_scale_is_ln   = sum(r$adj.P.Val < 0.05 & abs(r$logFC) >= cc * log(2)))))
}))
curve <- do.call(rbind, lapply(seq(0, 2, by = 0.05), function(cc) data.frame(
  log2FC_cutoff = cc,
  scale = c("read as log2 (factor 1)", "read as ln (factor 1.443)"),
  n = c(sum(res$as_is$adj.P.Val < 0.05 & abs(res$as_is$logFC) >= cc),
        sum(res$as_is$adj.P.Val < 0.05 & abs(res$as_is$logFC) >= cc * log(2))))))
write.csv(sens,  file.path(RESULTS, "sensitivity_counts.csv"), row.names = FALSE)
write.csv(curve, file.path(RESULTS, "sensitivity_curve.csv"),  row.names = FALSE)
write.csv(data.frame(check = names(inv), value = inv), file.path(RESULTS, "sensitivity_invariance.csv"), row.names = FALSE)
print(sens, row.names = FALSE); print(inv)
