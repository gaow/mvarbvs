# TO DO: Explain here what this script is for, and how to use it.
library(readr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(mvsusieR)
source("../code/mvsusie_plots.R")

# The ordering of the blood cell traits in the effect plots.  traits
blood_cell_traits <-
  c("RBC_count","Haemoglobin","MCV","RDW","MSCV","Reticulocyte_perc",
    "HLR_perc","Platelet_count", "Plateletcrit","PDW","WBC_count",
    "Lymphocyte_perc","Monocyte_perc","Neutrophill_perc","Eosinophill_perc",
    "Basophill_perc")

# Read the seq_gene data.
# seq_gene <- read_delim("../data/seq_gene.md.gz",delim = "\t",quote = "",
#                     col_types = cols(chromosome = "c"))
# class(seq_gene) <- "data.frame"
# seq_gene <- subset(seq_gene,
#                   group_label == "GRCh37.p5-Primary Assembly" &
#                    feature_type == "GENE")
# seq_gene <- transform(seq_gene,
#                       chr_start = chr_start/1e6,
#                       chr_stop = chr_stop/1e6)

# RUNX1 example.
poslim <- c(36.15,36.55)
dat <- readRDS(paste0("../output/blood_cell_traits/summary_stats/",
                      "bloodcells_chr21.36094353.36965761.summary_stats.rds"))
fit <- readRDS(paste0("../output/blood_cell_traits/mvsusie/",
                      "bloodcells_chr21.36094353.36965761.LDoriginal.Ycor.",
                      "mvsusierss.rds"))
p1 <- mvsusie_plot(fit,pos = dat$meta$POS/1e6,chr = 21,poslim = poslim,
                   conditions = blood_cell_traits)
p2 <- plot_gene_tracks(seq_gene,chr = 21,poslim = poslim,genes = "RUNX1")
#
# TO DO:
#
#   - Add SNP ids to effect plot.
#
#   - Reorder the traits (rows) in the effect plot.
#
#   - Vary color and size of dots in effect plot according to effects.
#
#   - Flip the effect directions when needed.
#
print(plot_grid(p1$pip_plot,p1$effect_plot,
                p2$plot,
                nrow = 2,ncol = 2,align = "v",axis = "lr",
                rel_heights = c(2,1),rel_widths = c(3,2)))
