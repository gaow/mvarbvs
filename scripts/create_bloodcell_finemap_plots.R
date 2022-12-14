# TO DO: Explain here what this script is for, and how to use it.
library(ggplot2)
library(ggrepel)
library(cowplot)
library(mvsusieR)
source("../code/mvsusie_plots.R")

# RUNX1 example.
region <- 
dat <- readRDS(paste0("../output/blood_cell_traits/summary_stats/",
                      "bloodcells_chr21.36094353.36965761.summary_stats.rds"))
fit <- readRDS(paste0("../output/blood_cell_traits/mvsusie/",
                      "bloodcells_chr21.36094353.36965761.LDoriginal.Ycor.",
                      "mvsusierss.rds"))
p1 <- pip_plot(fit,
               pos = dat$meta$POS/1e6,
               poslim = c(36.15,36.55))
#
# TO DO:
#
#  - Add ids of top ("sentinel") SNPs, and indicate the CS for these
#    sentinel SNPs.
#
#  - Add gene tracks.
# 
print(p1)
