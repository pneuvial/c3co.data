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

filenameSeg <- sprintf("segDat-%s.rds", patientID)
pathnameSeg <- file.path(path, filenameSeg)

if(!file.exists(pathnameSeg)){
  resC3co <- c3co(datList, parameters.grid = parameters.grid, verbose = TRUE)
  segdat <- resC3co@segDat
  segdat$bkp <- resC3co@bkp
  saveRDS(segdat, pathnameSeg)
}else{
  resC3co <- c3co(NULL, pathSeg = pathnameSeg, parameters.grid = parameters.grid, verbose = TRUE)
}
filenameResC3CO <- sprintf("resC3CO-%s.rds", patientID)
pathnameResC3CO <- file.path(path, filenameResC3CO)
saveRDS(resC3co, pathnameResC3CO)

