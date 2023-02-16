library(ggplot2)
library(cowplot)
library(ggrepel)
dat <- read.csv("../output/blood_cell_traits/cs1snp.csv",
                comment.char = "")
p <- ggplot(dat,aes(x = susie,y = mvsusie,label = trait,color = cell_type)) +
  geom_point(shape = 20) +
  geom_text_repel(size = 2.5,max.overlaps = Inf) +
  geom_abline(intercept = 0,slope = 1,color = "black",
              linetype = "dashed") +
  scale_x_continuous(limits = c(0,150),breaks = seq(0,400,50)) +
  scale_y_continuous(limits = c(0,335),breaks = seq(0,400,50)) +
  scale_color_manual(values = c("dodgerblue","salmon","red","violet")) +
  coord_fixed() +
  theme_cowplot(font_size = 12)
print(p)
ggsave("../plots/susie_vs_mvsusie_num_cs1snp.eps",p,height = 4,width = 4)
