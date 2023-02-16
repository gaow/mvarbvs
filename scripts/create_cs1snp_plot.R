library(ggplot2)
library(cowplot)
library(ggrepel)
dat <- read.csv("../output/blood_cell_traits/cs1snp.csv")
p <- ggplot(dat,aes(x = susie,y = mvsusie,label = trait)) +
  geom_point(shape = 21,color = "dodgerblue") +
  geom_text_repel(size = 2.5,color = "darkblue",max.overlaps = Inf) +
  geom_abline(intercept = 0,slope = 1,color = "black",
              linetype = "dashed") +
  scale_x_continuous(limits = c(0,150),breaks = seq(0,400,50)) +
  scale_y_continuous(limits = c(0,335),breaks = seq(0,400,50)) +
  coord_fixed() +
  theme_cowplot(font_size = 12)
