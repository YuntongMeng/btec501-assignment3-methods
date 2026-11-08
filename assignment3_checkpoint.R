# BTEC 501 Assignment 3 - Reproducible R checkpoint
# Run from the assignment directory that contains array_data.txt.
# In RStudio, run source("assignment3_checkpoint.R").
# From a terminal, run Rscript assignment3_checkpoint.R [input_file] [output_dir].
# This is the R walkthrough. The course handout specifies ExpressAnalyst
# for the submitted analysis.

run_assignment3 <- function(
    input_file = "array_data.txt",
    output_dir = "assignment3_r_results") {
  packages <- c("limma", "AnnotationDbi", "hgu133plus2.db", "org.Hs.eg.db",
                "clusterProfiler", "BiocParallel", "ggplot2", "pheatmap",
                "enrichplot")
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Install the missing R packages: ", paste(missing_packages, collapse = ", "))
  }
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file, "\nWorking directory: ", getwd())
  }

  # 1. Read the probe-set matrix and the group label for each sample.
  input <- read.delim(input_file, check.names = FALSE, colClasses = "character")
  if (nrow(input) < 3L || ncol(input) < 5L || names(input)[1] != "#NAME" ||
      input[1, 1] != "#CLASS") {
    stop("Expected #NAME in the header and #CLASS in the first data row.")
  }
  sample_ids <- names(input)[-1]
  probe_ids <- input[-1, 1]
  group <- factor(unlist(input[1, -1], use.names = FALSE), levels = c("1", "2"))
  expr <- as.matrix(input[-1, -1, drop = FALSE])
  storage.mode(expr) <- "double"
  dimnames(expr) <- list(probe_ids, sample_ids)
  if (anyNA(group) || any(table(group) < 2L) || any(!is.finite(expr)) ||
      anyNA(probe_ids) || any(!nzchar(probe_ids)) || anyDuplicated(probe_ids) ||
      any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop("Check for invalid groups, nonnumeric values, or missing/duplicate IDs.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output <- function(filename) file.path(output_dir, filename)
  save_plot <- function(filename, plot, width = 8, height = 5) {
    ggplot2::ggsave(output(filename), plot = plot, width = width,
                    height = height, dpi = 180, bg = "white")
  }

  fdr_cutoff <- 0.05
  signal_cutoffs <- c(0, 0.585, 1)
  strict_cutoff <- 1
  min_go_size <- 5L
  group_colors <- c("Group 1" = "#377EB8", "Group 2" = "#E64B35")
  group_label <- paste("Group", group)
  sample_table <- data.frame(sample_id = sample_ids, group = group_label)
  write.csv(sample_table, output("sample_groups.csv"), row.names = FALSE)

  # 2. Apply quantile normalization to the supplied expression matrix.
  # The course file does not document its earlier processing or log base.
  # The cutoffs below are differences on the supplied signal scale.
  normalized <- limma::normalizeBetweenArrays(expr, method = "quantile")
  qc_table <- rbind(
    data.frame(sample = rep(sample_ids, each = nrow(expr)),
               group = rep(group_label, each = nrow(expr)),
               signal = as.vector(expr), stage = "Before normalization"),
    data.frame(sample = rep(sample_ids, each = nrow(expr)),
               group = rep(group_label, each = nrow(expr)),
               signal = as.vector(normalized), stage = "After normalization")
  )
  qc_table$sample <- factor(qc_table$sample, levels = sample_ids)
  qc_table$stage <- factor(qc_table$stage,
                           levels = c("Before normalization", "After normalization"))
  boxplot <- ggplot2::ggplot(qc_table, ggplot2::aes(sample, signal, fill = group)) +
    ggplot2::geom_boxplot(outlier.size = 0.3) +
    ggplot2::facet_wrap(~stage, ncol = 1) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::labs(x = NULL, y = "Expression signal", fill = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
  save_plot("01_qc_boxplots.png", boxplot, width = 12, height = 7)

  density_plot <- ggplot2::ggplot(
    qc_table, ggplot2::aes(signal, group = sample, color = group)
  ) +
    ggplot2::geom_density(alpha = 0.45, linewidth = 0.4) +
    ggplot2::facet_wrap(~stage) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::labs(x = "Expression signal", y = "Density", color = NULL) +
    ggplot2::theme_bw()
  save_plot("02_qc_density.png", density_plot, width = 11)

  # The original PCA used the input matrix with centered, unscaled probe sets.
  pca <- prcomp(t(expr), center = TRUE, scale. = FALSE)
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pca_table <- data.frame(sample_table, PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  distances <- as.matrix(dist(pca_table[, c("PC1", "PC2")]))
  cross_group <- distances[group == "1", group == "2", drop = FALSE]
  closest <- which(cross_group == min(cross_group), arr.ind = TRUE)[1, ]
  closest_samples <- c(sample_ids[group == "1"][closest[1]],
                       sample_ids[group == "2"][closest[2]])
  axis_limits <- range(pca_table$PC1, pca_table$PC2)
  pca_plot <- ggplot2::ggplot(pca_table, ggplot2::aes(PC1, PC2, color = group)) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_text(
      data = pca_table[pca_table$sample_id %in% closest_samples, ],
      ggplot2::aes(label = sample_id), hjust = c(-0.15, 1.15), show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::coord_fixed(xlim = axis_limits, ylim = axis_limits) +
    ggplot2::labs(
      title = "Sample PCA before quantile normalization", color = NULL,
      x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2])
    ) +
    ggplot2::theme_bw()
  save_plot("03_qc_pca.png", pca_plot, width = 8.5, height = 7)
  write.csv(pca_table, output("pca_scores.csv"), row.names = FALSE)

  # 3. Fit the Group 2 minus Group 1 contrast and adjust across all probe sets.
  design <- model.matrix(~ group)
  fit <- limma::eBayes(limma::lmFit(normalized, design))
  de <- limma::topTable(fit, coef = "group2", number = Inf, adjust.method = "BH")
  de$probe_id <- rownames(de)
  names(de)[names(de) == "logFC"] <- "signal_difference"
  de <- de[, c("probe_id", setdiff(names(de), "probe_id"))]
  pass_fdr <- de$adj.P.Val < fdr_cutoff
  signal <- de$signal_difference
  threshold_counts <- do.call(rbind, lapply(signal_cutoffs, function(cutoff) {
    pass <- pass_fdr & abs(signal) >= cutoff
    data.frame(cutoff = cutoff, total = sum(pass),
               higher_in_group_2 = sum(pass & signal > 0),
               higher_in_group_1 = sum(pass & signal < 0))
  }))
  selected <- pass_fdr & abs(signal) >= strict_cutoff
  write.csv(de, output("differential_all_probes.csv"), row.names = FALSE)
  write.csv(de[selected, ], output("differential_strict_probes.csv"), row.names = FALSE)
  write.csv(threshold_counts, output("differential_threshold_counts.csv"), row.names = FALSE)

  for (cutoff in signal_cutoffs[-1]) {
    pass <- pass_fdr & abs(signal) >= cutoff
    plot_data <- data.frame(
      signal_difference = signal,
      minus_log10_fdr = -log10(pmax(de$adj.P.Val, .Machine$double.xmin)),
      result = ifelse(pass, ifelse(signal > 0, "Group 2", "Group 1"), "Other")
    )
    volcano <- ggplot2::ggplot(
      plot_data, ggplot2::aes(signal_difference, minus_log10_fdr, color = result)
    ) +
      ggplot2::geom_vline(xintercept = c(-cutoff, cutoff), linetype = "dashed") +
      ggplot2::geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed") +
      ggplot2::geom_point(size = 1.6, alpha = 0.8) +
      ggplot2::scale_color_manual(
        values = c(group_colors, Other = "grey75"),
        breaks = c("Group 1", "Group 2", "Other"),
        labels = c("Higher in Group 1", "Higher in Group 2", "Other")
      ) +
      ggplot2::labs(
        title = sprintf("Differential expression at signal cutoff %s", cutoff),
        x = "Signal difference (Group 2 - Group 1)",
        y = "-log10(BH-adjusted p)", color = NULL
      ) +
      ggplot2::theme_bw()
    save_plot(sprintf("04_volcano_cutoff_%s.png", cutoff), volcano)
  }

  # 4. Cluster the strictly selected probes with row Z scores.
  # The selected rows depend on the group comparison, so this is descriptive.
  selected_probes <- de$probe_id[selected]
  if (length(selected_probes) >= 2L) {
    heatmap_matrix <- t(scale(t(normalized[selected_probes, , drop = FALSE])))
    sample_tree <- hclust(dist(t(heatmap_matrix)), method = "complete")
    annotation <- data.frame(Group = group_label, row.names = sample_ids)
    pheatmap::pheatmap(
      heatmap_matrix, cluster_cols = sample_tree, clustering_method = "complete",
      annotation_col = annotation, annotation_colors = list(Group = group_colors),
      show_rownames = FALSE, border_color = NA, fontsize_col = 8,
      main = "Strictly selected probe sets (row Z scores)",
      filename = output("05_clustered_heatmap.png"), width = 11, height = 7
    )
    sample_table$cluster <- unname(cutree(sample_tree, k = 2)[sample_ids])
    write.csv(sample_table, output("heatmap_sample_clusters.csv"), row.names = FALSE)
  }

  # 5. Keep only unambiguous probe-to-Entrez mappings for both analyses.
  candidates <- AnnotationDbi::mapIds(
    hgu133plus2.db::hgu133plus2.db, keys = de$probe_id, keytype = "PROBEID",
    column = "ENTREZID", multiVals = "list"
  )
  candidates <- lapply(candidates[de$probe_id], function(ids) unique(ids[!is.na(ids)]))
  candidate_count <- lengths(candidates)
  entrez <- vapply(candidates, function(ids) {
    if (length(ids) == 1L) ids else NA_character_
  }, character(1))
  mapping <- data.frame(
    probe_id = de$probe_id, entrez_id = unname(entrez),
    candidate_count = unname(candidate_count),
    candidate_ids = vapply(candidates, paste, character(1), collapse = "/")
  )
  write.csv(mapping, output("probe_annotation.csv"), row.names = FALSE)
  background_genes <- unique(entrez[!is.na(entrez)])
  selected_genes <- unique(entrez[selected & !is.na(entrez)])
  org_db <- org.Hs.eg.db::org.Hs.eg.db

  # ORA uses the tested, uniquely annotated genes as the background.
  ora <- clusterProfiler::enrichGO(
    gene = selected_genes, universe = background_genes, OrgDb = org_db,
    keyType = "ENTREZID", ont = "BP", minGSSize = min_go_size,
    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1
  )
  ora_table <- as.data.frame(ora)
  if (nrow(ora_table)) {
    ora_table <- ora_table[order(ora_table$pvalue), , drop = FALSE]
  }
  ora_significant <- sum(ora_table$p.adjust < fdr_cutoff, na.rm = TRUE)
  write.csv(ora_table, output("go_bp_ora.csv"), row.names = FALSE)
  if (nrow(ora_table)) {
    top_ora <- head(ora_table, 6)
    top_ora$Description <- factor(top_ora$Description, levels = rev(top_ora$Description))
    ora_plot <- ggplot2::ggplot(
      top_ora, ggplot2::aes(FoldEnrichment, Description, size = Count)
    ) +
      ggplot2::geom_point(color = "grey35") +
      ggplot2::labs(
        title = "GO BP ORA - top terms by raw p-value",
        subtitle = sprintf("%d terms passed BH-adjusted p < %.2f", ora_significant, fdr_cutoff),
        x = "Fold enrichment", y = NULL, size = "Genes"
      ) +
      ggplot2::theme_bw()
    save_plot("06_ora_top_terms.png", ora_plot, width = 10)
  }

  # GSEA uses one signed moderated t per gene, without a differential cutoff.
  # Use the highest mean signal to choose a representative for duplicate genes.
  ranks <- data.frame(probe_id = de$probe_id, entrez_id = unname(entrez),
                      t_stat = de$t, mean_signal = de$AveExpr)
  ranks <- subset(ranks, !is.na(entrez_id) & is.finite(t_stat))
  ranks <- ranks[order(-ranks$mean_signal, ranks$probe_id), ]
  ranks <- ranks[!duplicated(ranks$entrez_id), ]
  ranks <- ranks[order(-ranks$t_stat, ranks$entrez_id), ]
  gene_rank <- setNames(ranks$t_stat, ranks$entrez_id)
  write.csv(ranks, output("gsea_gene_ranking.csv"), row.names = FALSE)

  previous_backend <- BiocParallel::bpparam()
  BiocParallel::register(BiocParallel::SerialParam())
  on.exit(BiocParallel::register(previous_backend), add = TRUE)
  set.seed(123)
  gsea <- clusterProfiler::gseGO(
    geneList = gene_rank, OrgDb = org_db, keyType = "ENTREZID",
    ont = "BP", minGSSize = min_go_size, pAdjustMethod = "BH",
    pvalueCutoff = 1, verbose = FALSE
  )
  gsea_table <- as.data.frame(gsea)
  if (nrow(gsea_table)) {
    gsea_table <- gsea_table[order(gsea_table$p.adjust, gsea_table$pvalue), , drop = FALSE]
  }
  gsea_significant <- sum(gsea_table$p.adjust < fdr_cutoff, na.rm = TRUE)
  write.csv(gsea_table, output("go_bp_gsea.csv"), row.names = FALSE)
  if (nrow(gsea_table)) {
    top_gsea <- gsea_table[1, ]
    gsea_plot <- enrichplot::gseaplot2(
      gsea, geneSetID = top_gsea$ID,
      title = sprintf("%s (BH-adjusted p = %.3g)", top_gsea$Description, top_gsea$p.adjust)
    )
    save_plot("07_gsea_top_term.png", gsea_plot, width = 9, height = 7)
  }

  summary <- c(
    sprintf("Input: %s", normalizePath(input_file, winslash = "/")),
    sprintf("%d probe sets; %d samples (Group 1: %d, Group 2: %d)",
            nrow(expr), ncol(expr), sum(group == "1"), sum(group == "2")),
    sprintf("Closest cross-group PCA pair: %s", paste(closest_samples, collapse = ", ")),
    "Contrast: Group 2 minus Group 1. Cutoffs are differences on the supplied signal scale.",
    capture.output(print(threshold_counts, row.names = FALSE)),
    sprintf("Annotation: %d unique, %d no match, %d ambiguous; %d distinct genes",
            sum(candidate_count == 1), sum(candidate_count == 0),
            sum(candidate_count > 1), length(background_genes)),
    sprintf("ORA: %d selected mapped genes, %d terms returned, %d pass BH p < %.2f",
            length(selected_genes), nrow(ora_table), ora_significant, fdr_cutoff),
    sprintf("GSEA: %d ranked genes, %d terms returned, %d pass BH p < %.2f",
            length(gene_rank), nrow(gsea_table), gsea_significant, fdr_cutoff),
    "GO BP only; minimum gene-set size 5. Interpret the enrichment within this teaching subset.",
    "", capture.output(sessionInfo())
  )
  writeLines(summary, output("checkpoint_summary.txt"))
  cat(paste(head(summary, 11), collapse = "\n"), "\n")
  message("Results saved to ", normalizePath(output_dir, winslash = "/"))
  invisible(list(differential = de, ora = ora, gsea = gsea, gene_rank = gene_rank))
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) > 2L) stop("Usage: Rscript assignment3_checkpoint.R [input_file] [output_dir]")
  do.call(run_assignment3, as.list(arguments))
} else {
  assignment3 <- run_assignment3()
}
