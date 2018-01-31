libraryBioc <- function(pkg) {
  if (require(pkg, character.only = TRUE)) return()
  source("https://bioconductor.org/biocLite.R")
  biocLite(pkg)
  library(pkg, character.only = TRUE)
}

R.utils::use("PSCBS")
libraryBioc("aroma.light")

## PSCBS data on GSE22298 perform PSCBS data
path <- file.path("locusData", "GSE22928,MM")
pathSegPSCBS <- "PSCBSdata"

pathnames <- list.files(path, pattern = "[.]rds$", full.names = TRUE)
datToSeg <- lapply(pathnames, FUN = function(ff) {
  df <- readRDS(ff)
  ## Chromosomes 1 to 2
  subset(df, chromosome %in% 1:2)
})
names(datToSeg) <- tools::file_path_sans_ext(basename(pathnames))


### Run PSCBS (5 samples and 2 chromosomes)
PSCBSdata <- lapply(datToSeg[1:5], FUN = function(df) {
  gaps <- findLargeGaps(df, minLength = 1e+06)
  knownSegments <- gapsToSegments(gaps)
  segmentByPairedPSCBS(df, knownSegments = knownSegments, preserveScale = FALSE, seed = 48879, verbose = -10)
})
for (name in names(PSCBSdata)) sampleName(PSCBSdata[[name]]) <- name

stopifnot(length(PSCBSdata) == 5)
stopifnot(!is.null(names(PSCBSdata)))
saveRDS(PSCBSdata, file = sprintf("%s.rds", "GSE22928,MM,1-5"), compress = "bzip2")

## Save
devtools::use_data(PSCBSdata)
