##########################################################################
# Run c3co for patient RK29 from GSE47077
##########################################################################

library("c3co")
patientID <- "RK29"
stopifnot(packageVersion("c3co")>='0.2.1.9000')

path <- "GSE47077"
path <- R.utils::Arguments$getReadablePath(path)
filename <- sprintf("dat-%s.rds", patientID)
pathname <- file.path(path, filename)
datList <- readRDS(pathname)

lambda.grid <- 10^(-seq(from = 2, to = 6, by = 1)) ## penalty
p.list <- 2:length(datList) ## candidate number of subclones
parameters.grid <- list(lambda = lambda.grid, nb.arch = p.list)

filenameSeg <- sprintf("segDat.rda", patientID)
pathnameSeg <- file.path(path, filenameSeg)

if(!file.exists(pathnameSeg)){
  resC3CO <- c3co(datList, parameters.grid = parameters.grid, verbose = TRUE)
  segDat <- resC3CO@segDat
  segDat$bkp <- resC3CO@bkp
  devtools::use_data(segDat, path)
}else{
  resC3CO <- c3co(NULL, pathSeg = pathnameSeg, parameters.grid = parameters.grid, verbose = TRUE)
}
devtools::use_data(resC3CO, path)

