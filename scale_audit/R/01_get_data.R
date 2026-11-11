# 01  Download the two public inputs:
#   (a) MAMA 2.2.1 source from the CRAN archive -> ColonData, "Japan" component
#   (b) GEO GSE4554 series matrix -> original MAS5 signals for all 84 arrays
source("scale_audit/R/common.R")
suppressPackageStartupMessages(library(GEOquery))

# (a) MAMA ColonData ----------------------------------------------------------
tgz <- file.path(DATA, "MAMA_2.2.1.tar.gz")
if (!file.exists(tgz))
  download.file("https://cran.r-project.org/src/contrib/Archive/MAMA/MAMA_2.2.1.tar.gz", tgz, mode = "wb")
untar(tgz, files = "MAMA/data/ColonData.RData", exdir = DATA)
env <- new.env()
load(file.path(DATA, "MAMA/data/ColonData.RData"), envir = env)
cd <- env$ColonData                       # S4 "MetaArray"; read slots as attributes
names_ds <- attr(cd, "datanames")
japan    <- attr(cd, "GEDM")[[which(names_ds == "Japan")]]
japan_cl <- attr(cd, "clinical")[[which(names_ds == "Japan")]]
rownames(japan_cl) <- colnames(japan)
other    <- setNames(attr(cd, "GEDM")[names_ds != "Japan"], names_ds[names_ds != "Japan"])

# (b) GEO GSE4554 -------------------------------------------------------------
es   <- getGEO("GSE4554", GSEMatrix = TRUE, destdir = DATA)[[1]]
mas5 <- exprs(es)                         # MAS5 signal, not log-transformed
pd   <- pData(es)

saveRDS(list(japan = japan, japan_cl = japan_cl, other = other,
             mas5 = mas5, pdata = pd), file.path(DATA, "inputs.rds"))
cat("MAMA Japan:", dim(japan), "| GSE4554 MAS5:", dim(mas5), "\n")
