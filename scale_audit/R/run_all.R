# Run the whole audit from the repository root:  Rscript scale_audit/R/run_all.R
for (s in c("01_get_data", "02_provenance", "03_scale_tests", "04_sensitivity", "05_figures")) {
  message("== ", s)
  source(file.path("scale_audit/R", paste0(s, ".R")), local = new.env())
}
writeLines(capture.output(sessionInfo()), "scale_audit/results/sessionInfo.txt")
