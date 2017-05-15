## PSCBS data on GSE22298 perform PSCBS data
path <- "locusData/GSE22928,MM/"
pathSegPSCBS <- "PSCBSdata"

datToSeg <- lapply(list.files(path), function (ff) {
  df <- readRDS(file.path(path, ff))
  ## Chromosomes 1 to 2
  df <- subset(df, chromosome%in%1:2)
})

## Name of each sample
names(datToSeg) <- list.files(path)

install.packages("PSCBS")
source("https://bioconductor.org/biocLite.R")
biocLite("aroma.light")

### Run PSCBS (test on 5 samples and 5 chromosomes)
PSCBSdata <- lapply(1:5, function (ff) {
  df <- datToSeg[[ff]]
  gaps <- PSCBS::findLargeGaps(df, minLength = 1e+06)
  knownSegments <- PSCBS::gapsToSegments(gaps)
  fit <- PSCBS::segmentByPairedPSCBS(df, knownSegments = knownSegments, preserveScale = FALSE, seed = 48879, verbose = -10)
  return(fit)
})
## Save
devtools::use_data(PSCBSdata)

