# TO DO: Explain here what this script is for, and how to use it.
library(susieR)
library(mvsusieR)
susie <- readRDS("../output/blood_cell_traits/susie/EXT1_susie.rds")
mvsusie <- readRDS("../output/blood_cell_traits/mvsusie/EXT1_mvsusie.rds")
print(mvsusie$sets$purity)
j <- mvsusie$sets$cs$L1
print(mvsusie$pip[j])
print(which(sapply(susie,function (x) !is.null(x$sets$cs))))
