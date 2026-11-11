# Shared paths and helpers for the scale audit.
# Run every script from the repository root, e.g.  Rscript scale_audit/R/run_all.R

suppressPackageStartupMessages({
  library(Biobase)
})

AUDIT   <- "scale_audit"
DATA    <- file.path(AUDIT, "data")      # downloaded inputs (git-ignored)
RESULTS <- file.path(AUDIT, "results")   # small CSV tables (committed)
FIGS    <- file.path(AUDIT, "figures")   # PNG figures (committed)
for (d in c(DATA, RESULTS, FIGS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Categorical colours (validated: lightness, chroma, CVD separation)
COL <- c(ln = "#2a78d6", log2 = "#eb6834", glog2 = "#1baf7a")
INK <- "#3b3b38"; GRID <- "#e6e5e0"

theme_audit <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = GRID, linewidth = 0.3),
                   text = ggplot2::element_text(colour = INK),
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                   plot.title.position = "plot",
                   legend.position = "top")
}

# Plain quantile normalisation via preprocessCore (as named in the MAMA author's thesis)
qn <- function(M) {
  out <- preprocessCore::normalize.quantiles(as.matrix(M))
  dimnames(out) <- dimnames(M); out
}

load_inputs <- function() {
  f <- file.path(DATA, "inputs.rds")
  if (!file.exists(f)) stop("Run scale_audit/R/01_get_data.R first.")
  readRDS(f)
}
