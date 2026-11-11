# 05  Figures for the report (PNG, 300 dpi)
source("scale_audit/R/common.R")
suppressPackageStartupMessages(library(ggplot2))
x    <- load_inputs()
J    <- x$japan
cand <- readRDS(file.path(DATA, "candidates.rds"))
fits <- read.csv(file.path(RESULTS, "scale_fits_per_sample.csv"))
summ <- read.csv(file.path(RESULTS, "scale_summary.csv"))
curve <- read.csv(file.path(RESULTS, "sensitivity_curve.csv"))
lab  <- c(ln = "H1: ln + QN", log2 = "H2: log2 + QN", glog2 = "H3: VSN (glog2) + QN")
save <- function(p, f, w, h) ggsave(file.path(FIGS, f), p, width = w, height = h, dpi = 300, bg = "white")

# Figure 1  provenance ---------------------------------------------------------
boxes <- data.frame(
  x = 1:4,
  head = c("GEO GSE4554", "Thesis, section 5.2", "MAMA ColonData", "This audit"),
  body = c("84 colorectal tumours\nHG-U133 Plus 2.0\nMAS5 signal (not logged)",
           "MAS5 -> VSN ->\nnormalize.quantiles\n(Ihnatova, 2010)",
           "\"Japan\": 41 samples\n500 probe sets\ndocumented as log2",
           "compare with MAS5\nthree scale hypotheses\nsensitivity of DE"))
p1 <- ggplot(boxes) +
  geom_segment(data = data.frame(x = 1:3 + 0.36, xend = 2:4 - 0.36), aes(x = x, xend = xend, y = 0, yend = 0),
               arrow = arrow(length = unit(2.5, "mm")), colour = INK, linewidth = 0.5) +
  geom_label(aes(x = x, y = 0.28, label = head), fontface = "bold", size = 3.4, label.size = 0, colour = INK, fill = "white") +
  geom_label(aes(x = x, y = -0.05, label = body), size = 3, lineheight = 1.05, label.size = 0.3,
             label.padding = unit(2.5, "mm"), colour = INK, fill = "#f4f3ef") +
  scale_x_continuous(limits = c(0.5, 4.5)) + scale_y_continuous(limits = c(-0.22, 0.36)) +
  theme_void()
save(p1, "fig1_provenance.png", 9, 1.7)

# Figure 2  shape vs value for one sample ---------------------------------------
s <- "GSM101849"
d2 <- do.call(rbind, lapply(names(cand), function(h)
  data.frame(hypothesis = factor(lab[h], lab), candidate = cand[[h]][, s], mama = J[, s], row.names = NULL)))
ann <- do.call(rbind, lapply(names(cand), function(h) {
  f <- fits[fits$hypothesis == h & fits$GSM == s, ]
  data.frame(hypothesis = factor(lab[h], lab), txt = sprintf("R-squared = %.3f\nslope = %.2f", f$R2, f$slope))
}))
p2 <- ggplot(d2, aes(candidate, mama)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#8a8983", linewidth = 0.5) +
  geom_point(aes(colour = hypothesis), size = 0.9, alpha = 0.55, show.legend = FALSE) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = INK, linewidth = 0.6) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = txt), hjust = -0.1, vjust = 1.2, size = 3.1, colour = INK) +
  facet_wrap(~hypothesis, nrow = 1) +
  scale_colour_manual(values = setNames(COL, lab)) +
  coord_equal(xlim = c(-1.5, 14), ylim = c(-1.5, 14)) +
  labs(x = "Value predicted by hypothesis", y = "MAMA Japan value",
       title = paste("Shape vs. value, sample", s),
       subtitle = "Dashed line: identical values (y = x). Solid line: free linear fit. All three fit the shape; only H1 lies near y = x.") +
  theme_audit()
save(p2, "fig2_shape_vs_value.png", 9, 3.8)

# Figure 3  per-sample slopes --------------------------------------------------
fits$hyp <- factor(lab[fits$hypothesis], rev(lab))
p3 <- ggplot(fits, aes(slope, hyp, colour = hypothesis)) +
  geom_vline(xintercept = 1, linetype = "22", colour = "#8a8983") +
  geom_jitter(height = 0.12, width = 0, size = 1.8, alpha = 0.8, show.legend = FALSE) +
  geom_text(data = transform(summ, hyp = factor(lab[hypothesis], rev(lab))),
            aes(x = 1.06, y = hyp, label = sprintf("median R-squared %.4f", R2_median)),
            hjust = 0, colour = INK, size = 3.1, inherit.aes = FALSE) +
  scale_colour_manual(values = COL) +
  scale_x_continuous(limits = c(0.55, 1.36), breaks = seq(0.6, 1.2, 0.1)) +
  labs(x = "Per-sample slope of MAMA value on hypothesis (1 = same scale)", y = NULL,
       title = "Scale agreement across the 41 samples",
       subtitle = "Identical shape fit (R-squared) for H1 and H2: a change of log base is a pure rescaling.") +
  theme_audit()
save(p3, "fig3_slopes.png", 8, 3.2)

# Figure 4  sensitivity of DE counts ------------------------------------------
curve$scale <- factor(curve$scale, c("read as log2 (factor 1)", "read as ln (factor 1.443)"))
ends <- do.call(rbind, lapply(split(curve, curve$scale), function(d) d[d$log2FC_cutoff == 1, ]))
p4 <- ggplot(curve, aes(log2FC_cutoff, n, colour = scale)) +
  geom_vline(xintercept = c(0.585, 1), colour = GRID, linewidth = 0.6) +
  geom_line(linewidth = 0.8) +
  geom_point(data = ends, size = 2.6) +
  geom_text(data = ends, aes(label = n), nudge_x = 0.07, hjust = 0, size = 3.3, colour = INK, show.legend = FALSE) +
  scale_colour_manual(values = unname(COL[c("log2", "ln")]), name = NULL) +
  labs(x = "|log2 fold change| cutoff", y = "Probe sets with adj. P < 0.05",
       title = "Differential expression (MSS vs MSI) under the two readings",
       subtitle = "Adjusted P-values are identical; only the fold-change cutoff moves.") +
  theme_audit()
save(p4, "fig4_sensitivity.png", 7, 4)
cat("Figures written to", FIGS, "\n")
