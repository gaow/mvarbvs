# Short script to inspect the results in more detail for the
# fine-mapping examples included in the main text.
library(susieR)
library(mvsusieR)
susie <- readRDS("../output/blood_cell_traits/susie/EXT1_susie.rds")
mvsusie <- readRDS("../output/blood_cell_traits/mvsusie/EXT1_mvsusie.rds")
print(mvsusie$sets$purity)
j <- mvsusie$sets$cs$L1
print(mvsusie$pip[j])
print(which(sapply(susie,function (x) !is.null(x$sets$cs))))
print(length(susie$Basophill_perc$sets$cs$L1))
