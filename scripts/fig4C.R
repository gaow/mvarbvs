library(dplyr)
library(ggplot2)
library(cowplot)

load("../output/blood_cell_traits/Fig4data.RData")

trait <-
  c("WBC_count", "RBC_count","Haemoglobin","MCV","RDW","Platelet_count",
    "Plateletcrit","PDW","Lymphocyte_perc","Monocyte_perc","Neutrophill_perc",
    "Eosinophill_perc","Basophill_perc","Reticulocyte_perc","MSCV","HLR_perc")
bloodcells_col <-
  cbind(trait,
        c("Compound white cell", "Mature red cell", "Mature red cell",
          "Mature red cell", "Mature red cell", "Platelet", "Platelet",
          "Platelet", "Compound white cell", "Compound white cell",
          "Compound white cell", "Compound white cell", "Compound white cell",
          "Immature red cell", "Mature red cell","Immature red cell"),
        c("#33cccc","red","red","red","red",
          "#cc66ff","#cc66ff","#cc66ff",
          "#33cccc","#33cccc","#33cccc","#33cccc","#33cccc",
          "pink","red","pink"))

rename <- list("WBC_count" = "WBC#",
              "RBC_count" = "RBC#",
              "Haemoglobin" = "HGB",
              "MCV" = "MCV",
              "RDW" = "RDW",
              "Platelet_count" = "PLT#",
              "Plateletcrit" = "PCT",
              "PDW" = "PDW",
              "Lymphocyte_perc" = "LYMPH%",
              "Monocyte_perc" = "MONO%",
              "Neutrophill_perc" = "NEUT%",
              "Eosinophill_perc" = "EO%",
              "Basophill_perc" = "BASO%",
              "Reticulocyte_perc" = "RET%",
              "MSCV" = "MSCV",
              "HLR_perc" = "HLR%")

CS_compare <-
  res_rss_region %>%
  group_by(trait) %>%
  summarize(CSnumber = sum(CS_num))
tmp <- colSums(res_mvrss_CS[,6:21])
CS_compare$mvCSnumber <- tmp[match(CS_compare$trait,names(tmp))]
CS_compare$type <- bloodcells_col[match(CS_compare$trait,bloodcells_col[,1]),2]
CS_compare$trait <- sapply(CS_compare$trait,function(x) rename[[x]])
class(CS_compare) <- "data.frame"
p <- ggplot(CS_compare,aes(x = CSnumber,y = mvCSnumber,color = type)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("#33cccc","pink","red","#cc66ff")) +
  geom_abline(slope = 1,intercept = 0,linetype = 2) +
    geom_text(aes(label = trait),vjust = "inward",hjust = "inward",size = 2) +
    labs(x = "SuSiE", y = "mvSuSiE") +
    scale_x_continuous(limits = c(125,800),breaks = seq(0,2000,200)) +
    scale_y_continuous(limits = c(125,1425),breaks = seq(0,2000,200)) +
    coord_fixed() +
    theme_cowplot(font_size = 12) +
    theme(legend.position = "none",
          title = element_text(size = 15))
print(p)
ggsave("susie_vs_mvsusie_num_cs.eps",p,height = 4,width = 4)
