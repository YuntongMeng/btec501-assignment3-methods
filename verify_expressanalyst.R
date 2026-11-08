# Reproduce the ExpressAnalyst differential expression comparison in R.
# Rscript verify_expressanalyst.R <expression.tsv> <gene_mapping_report.csv> [output.csv]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Usage: Rscript verify_expressanalyst.R <expression.tsv> <gene_mapping_report.csv> [output.csv]")
}

expression <- read.delim(args[1], check.names = FALSE, colClasses = "character")
if (nrow(expression) < 3L || ncol(expression) < 5L ||
    names(expression)[1] != "#NAME" || expression[1, 1] != "#CLASS") {
  stop("Expected a #NAME header and a #CLASS row with sample group labels.")
}

probe_ids <- expression[-1, 1]
sample_ids <- names(expression)[-1]
group <- factor(unlist(expression[1, -1], use.names = FALSE), levels = c("1", "2"))
signals <- as.matrix(expression[-1, -1, drop = FALSE])
storage.mode(signals) <- "double"
dimnames(signals) <- list(probe_ids, sample_ids)
if (anyNA(group) || any(table(group) < 2L) || any(!is.finite(signals)) ||
    anyNA(probe_ids) || any(!nzchar(probe_ids)) || anyDuplicated(probe_ids) ||
    any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
  stop("Check the group labels, sample and probe IDs, and expression values.")
}

mapping <- read.csv(args[2], colClasses = "character")
required_columns <- c("original_id", "entrez_id", "match_status")
if (!all(required_columns %in% names(mapping)) ||
    anyNA(mapping$original_id) || anyDuplicated(mapping$original_id) ||
    !setequal(mapping$original_id, probe_ids)) {
  stop("The mapping report must contain one row for each input probe set.")
}
mapping <- mapping[match(probe_ids, mapping$original_id), ]
if (any(!mapping$match_status %in% c("TRUE", "FALSE"))) {
  stop("Expected TRUE or FALSE in the mapping report's match_status column.")
}
matched <- mapping$match_status == "TRUE"
if (anyNA(mapping$entrez_id[matched]) || any(!grepl("^[0-9]+$", mapping$entrez_id[matched]))) {
  stop("Matched probes must have numeric Entrez IDs in the mapping report.")
}

# Average probes that map to the same gene within each sample. Retain
# unmatched probes under their original IDs, matching the website settings.
feature_ids <- ifelse(matched, mapping$entrez_id, probe_ids)
features <- limma::avereps(signals, ID = feature_ids)
normalized <- preprocessCore::normalize.quantiles(features)
dimnames(normalized) <- dimnames(features)

fit <- limma::eBayes(limma::lmFit(normalized, model.matrix(~ group)))
de <- limma::topTable(fit, coef = "group2", number = Inf, adjust.method = "BH")

output_file <- if (length(args) == 3L) args[3] else file.path("output", "expressanalyst_de.csv")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(GeneID = rownames(de), de), output_file, row.names = FALSE)

cutoffs <- c(0, 0.58, 1)
pass_fdr <- de$adj.P.Val < 0.05
counts <- data.frame(
  cutoff = cutoffs,
  total = vapply(cutoffs, function(x) sum(pass_fdr & abs(de$logFC) >= x), integer(1)),
  higher_in_group_2 = vapply(cutoffs, function(x) sum(pass_fdr & de$logFC > 0 & de$logFC >= x), integer(1)),
  higher_in_group_1 = vapply(cutoffs, function(x) sum(pass_fdr & de$logFC < 0 & de$logFC <= -x), integer(1))
)
print(counts, row.names = FALSE)
