libraryBioc <- function(pkg) {
  if (require(pkg, character.only = TRUE)) return()
  source("https://bioconductor.org/biocLite.R")
  biocLite(pkg)
  library(pkg, character.only = TRUE)
}

R.utils::use("R.utils")
use("R.filesets")

libraryBioc("Biobase")
libraryBioc("GEOquery")

## load series and platform data from GEO
dataSet <- "GSE22928"
filePath <- file.path(dataSet, list.files(dataSet, pattern="matrix_processed"))[1]
if(!file.exists(filePath)){
  message("load GSE22928 data set from GEO. Could take time")
  gset <- getGEOSuppFiles('GSE22928')
  filePath <- file.path(dataSet, list.files(dataSet, pattern="matrix_processed"))[1]
}
dataSNP <- readDataFrame(filePath)[, -(92:97)]

platform <- "GPL6984"
message(sprintf("load %s data set from GEO. Could take time", platform))
gsetGPL <- getGEO(platform, destdir=".")
posData <- gsetGPL@dataTable@table

message(sprintf("Merge the two data frame posData and data SNP"))
data <- merge(dataSNP, posData, by.x="ID_REF", by.y="ID")
data$Chr <- as.numeric(data$Chr)
o <- order(data$Chr, data$MapInfo)
dataOrdered <- data[o, ]

## annotation data: chromosome, position
chrC <- data[["Chr"]]
chr <- chrC
print(table(chr))


## order by genomic position
pos <- data[["MapInfo"]]
annDat <- cbind(chromosome=chr, x=pos)
print(head(annDat))

o <- order(annDat[, "chromosome"], annDat[, "x"])
annDat <- annDat[o, ]
print(head(annDat))

data <- data[o, ]
print(head(data))
rm(o)

n <- nrow(data)


## LRR and BAF
sampleNames <- gsub("_", "", colnames(data)[seq(from=2, to=91, by=6)])
## export paths
lpath <- "locusData"
lpath <- Arguments$getWritablePath(lpath)
## sample-specific paths
ds <- sprintf("%s,%s", dataSet, "MM")

slpath <- file.path(lpath, ds)
slpath <- Arguments$getWritablePath(slpath)
## - - - - - - - - - - - - - - - -
## 1- Save locusData
## - - - - - - - - - - - - - - - -

for (ss in seq_along(sampleNames)) {
  sampleName <- sampleNames[ss]
  message(sprintf("- Sample %d ('%s')", ss, sampleName))

  filename <- sprintf("%s,%s.rds", ds, sampleName)
  pathname <- file.path(slpath, filename)

  ## Already processed?
  if (file_test("-f", pathname)) next
  
  colsLR <- grep("Log_R_Ratio", colnames(data))[ss]
  colsBAF <- grep("B_Allele_Freq", colnames(data))[ss]
  y <- data[, c(colsLR, colsBAF)]
  idxBB <- which(y[, 2] > 0.9)
  idxAA <- which(y[, 2] < 0.1)
  idxAB <- which(y[, 2] >= 0.1 & y[, 2] <= 0.9)
  geno <- numeric(length=length(y[, 2]))
  geno[idxBB] <- 1
  geno[idxAA] <- 0
  geno[idxAB] <- 1/2
  
  betaT <- y[, 2]
  ## True BAF corresponding to genotypes
  gammaN <- geno
  betaN <- gammaN
  datN <- cbind(betaN, gammaN)
  names(y) <- c("CT", "betaT")
  y$CT <- 2*2^y$CT
  ## tumorBoost
  datPP <- cbind(annDat, y, betaN=betaN, betaTN=betaT, gammaN=gammaN)
  ## write to disk
  saveRDS(datPP, file=pathname)
}
