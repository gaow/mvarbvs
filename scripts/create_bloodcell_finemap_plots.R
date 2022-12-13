# TO DO: Explain here what this script is for, and how to use it.
library(ggplot2)
library(ggrepel)
library(cowplot)
library(mvsusieR)

# Colors from colorbrewer2.org.
cs_colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c",
               "#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928")

# RUNX1 example.
region <- c(36.15,36.55)
dat <- readRDS(paste0("../output/blood_cell_traits/summary_stats/",
                      "bloodcells_chr21.36094353.36965761.summary_stats.rds"))
fit <- readRDS(paste0("../output/blood_cell_traits/mvsusie/",
                      "bloodcells_chr21.36094353.36965761.LDoriginal.Ycor.",
                      "mvsusierss.rds"))
pdat1 <- data.frame(pip = fit$pip,
                    pos = dat$meta$POS/1e6,
                    cs  = 0)
L <- length(fit$sets$cs)
for (i in 1:L) {
  j <- fit$sets$cs[[i]]
  pdat1[j,"cs"] <- i
}
pdat2 <- subset(pdat1,cs > 0)
pdat2 <- transform(pdat2,cs = factor(cs))
pdat1 <- subset(pdat1,pos >= region[1] & pos <= region[2])
pdat2 <- subset(pdat2,pos >= region[1] & pos <= region[2])
p1 <- ggplot(pdat1,aes(x = pos,y = pip)) +
  geom_point(color = "darkblue",shape = 20,size = 1.25) +
  geom_point(shape = 1,size = 1.25,stroke = 1.25,
             data = pdat2,
             mapping = aes(x = pos,y = pip,color = cs)) +
  scale_color_manual(values = cs_colors) +
  guides(colour = guide_legend(override.aes = list(shape = 20,size = 1.5))) +
  labs(x = "chromosome 21 position (Mb)",y = "PIP",color = "CS") +
  theme_cowplot(font_size = 10)

# TO DO:
#
#  - Add ids of top SNPs, and indicate the CS.
#
#  - Add purity values, CS sizes.
#
#  - Add gene tracks.
# 
