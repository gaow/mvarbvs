# Script to generate the plots for the main fine-mapping examples in
# the paper.
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
seq_gene <- read_delim("../data/seq_gene.md.gz",delim = "\t",quote = "",
                      col_types = cols(chromosome = "c"))
class(seq_gene) <- "data.frame"
seq_gene <- subset(seq_gene,
                   group_label == "GRCh37.p5-Primary Assembly" &
                   feature_type == "GENE")
seq_gene <- transform(seq_gene,
                      chr_start = chr_start/1e6,
                      chr_stop = chr_stop/1e6)

# EXT1/SAMD12 example.
poslim <- c(118.81,119.675)
dat <- readRDS(paste0("../output/blood_cell_traits/summary_stats/",
                      "bloodcells_chr8.118878055.119378055.summary_stats.rds"))
fit <- readRDS(paste0("../output/blood_cell_traits/mvsusie/",
                      "bloodcells_chr8.118878055.119378055.LDoriginal.Ycor.",
                      "mvsusierss.rds"))
p1 <- mvsusie_plot(fit,pos = dat$meta$POS/1e6,markers = dat$meta$ID,chr = 8,
                   poslim = poslim,conditions = blood_cell_traits)
p2 <- plot_gene_tracks(seq_gene,chr = 8,poslim = poslim,
                       genes = c("EXT1","SAMD12"))
print(plot_grid(p1$pip_plot,p1$effect_plot,
                p2$plot,
                nrow = 2,ncol = 2,align = "v",axis = "lr",
                rel_heights = c(3,1),rel_widths = c(2,1)))
ggsave("../plots/bloodcells_finemap_ext1_samd12_pips.eps",
       plot_grid(p1$pip_plot,p2$plot,nrow = 2,ncol = 1,align = "v",axis = "lr",
                 rel_heights = c(2,1)),
       height = 2.75,width = 7)
ggsave("../plots/bloodcells_finemap_ext1_samd12_effects.pdf",
       p1$effect_plot,height = 3.25,width = 3)

# TNS3 example.
poslim <- c(47.1,47.7)
dat <- readRDS(paste0("../output/blood_cell_traits/summary_stats/",
                      "bloodcells_chr7.47187274.47687274.summary_stats.rds"))
fit <- readRDS(paste0("../output/blood_cell_traits/mvsusie/",
                      "bloodcells_chr7.47187274.47687274.LDoriginal.Ycor.",
                      "mvsusierss.rds"))
p1 <- mvsusie_plot(fit,pos = dat$meta$POS/1e6,markers = dat$meta$ID,chr = 7,
                   poslim = poslim,conditions = blood_cell_traits)
p2 <- plot_gene_tracks(seq_gene,chr = 7,poslim = poslim,genes = "TNS3")
print(plot_grid(p1$pip_plot,p1$effect_plot,
                p2$plot,
                nrow = 2,ncol = 2,align = "v",axis = "lr",
                rel_heights = c(3,1),rel_widths = c(2,1)))
ggsave("../plots/bloodcells_finemap_tns3_pips.eps",
       plot_grid(p1$pip_plot,p2$plot,nrow = 2,ncol = 1,align = "v",axis = "lr",
                 rel_heights = c(2,1)),
       height = 2.75,width = 7)
ggsave("../plots/bloodcells_finemap_tsn3_effects.pdf",
       p1$effect_plot,height = 3.25,width = 2.75)

# RUNX1 example.
poslim <- c(36.15,36.55)
dat <- readRDS(paste0("../output/blood_cell_traits/summary_stats/",
                      "bloodcells_chr21.36094353.36965761.summary_stats.rds"))
fit <- readRDS(paste0("../output/blood_cell_traits/mvsusie/",
                      "bloodcells_chr21.36094353.36965761.LDoriginal.Ycor.",
                      "mvsusierss.rds"))
p1 <- mvsusie_plot(fit,pos = dat$meta$POS/1e6,markers = dat$meta$ID,chr = 21,
                   poslim = poslim,conditions = blood_cell_traits)
p2 <- plot_gene_tracks(seq_gene,chr = 21,poslim = poslim,genes = "RUNX1")
print(plot_grid(p1$pip_plot,p1$effect_plot,
                p2$plot,
                nrow = 2,ncol = 2,align = "v",axis = "lr",
                rel_heights = c(3,1),rel_widths = c(3,2)))
ggsave("../plots/bloodcells_finemap_runx1_pips.eps",
       plot_grid(p1$pip_plot,p2$plot,nrow = 2,ncol = 1,align = "v",axis = "lr",
                 rel_heights = c(2,1)),
       height = 3.25,width = 7)
ggsave("../plots/bloodcells_finemap_runx1_effects.pdf",
       p1$effect_plot,height = 3.25,width = 4)

# ZFPM1-PIEZO1 example.
poslim <- c(88.45,89)
dat <- readRDS(paste0("../output/blood_cell_traits/summary_stats/",
                      "bloodcells_chr16.88306007.90291983.summary_stats.rds"))
fit <- readRDS(paste0("../output/blood_cell_traits/mvsusie/",
                      "bloodcells_chr16.88306007.90291983.LDoriginal.Ycor.",
                      "mvsusierss.rds"))
p1 <- mvsusie_plot(fit,pos = dat$meta$POS/1e6,markers = dat$meta$ID,chr = 16,
                   poslim = poslim,conditions = blood_cell_traits)
p2 <- plot_gene_tracks(seq_gene,chr = 16,poslim = poslim,
                       genes = c("ZNF469","ZFPM1","ZC3H18","IL17C","CYBA",
                                 "MVD","SNAI3","RNF166","CTU2",
                                 "PIEZO1","CDT1","APRT","GALNS",
                                 "TRAPPC2L","PABPN1L","CBFA2T3"))
print(plot_grid(p1$pip_plot,p1$effect_plot,
                p2$plot,
                nrow = 2,ncol = 2,align = "v",axis = "lr",
                rel_heights = c(3,1),rel_widths = c(3,2)))
ggsave("../plots/bloodcells_finemap_zfpm1_piezo1_pips.eps",
       plot_grid(p1$pip_plot,p2$plot,nrow = 2,ncol = 1,align = "v",axis = "lr",
                 rel_heights = c(2,1)),
       height = 3.25,width = 7)
ggsave("../plots/bloodcells_finemap_zfpm1_piezo1_effects.pdf",
       p1$effect_plot,height = 3.25,width = 4)
