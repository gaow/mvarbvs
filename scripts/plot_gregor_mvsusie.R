# TO DO: Explain here what this script is for, and how to use it.
library(ggplot2)
library(cowplot)
traits <- c("crosstrait",
            "RBC","HGB","MCV","RDW","MSCV",          # mature red blood cell
            "Reticulocyte_perc","HLR_perc",          # immature red blood cell
            "Platelet_count","Plateletcrit","PDW",   # platelet
            "WBC","Lymphocyte_perc","Monocyte_perc", # white blood cell
            "Neutrophill_perc","Eosinophill_perc","Basophill_perc")

# Repeat for each trait.
n <- length(traits)
scatterplots <- vector("list",n)
names(scatterplots) <- traits
for (trait in traits) {

  # Load the susie and mvsusie results.
  susie_file <- paste("../output/blood_cell_traits/gregor/SuSiE",trait,
                      "positive_cs.index_gregor_output_StatisticSummaryFile",
                      "fisher.txt",sep = "_")
  mvsusie_file <- paste("../output/blood_cell_traits/gregor/mvSuSiE",trait,
                        "positive_cs.index_gregor_output_StatisticSummaryFile",
                        "fisher.txt",sep = "_")
  susie   <- read.table(susie_file,sep = " ",header = TRUE,
                        stringsAsFactors = FALSE)
  mvsusie <- read.table(mvsusie_file,sep = " ",header = TRUE,
                        stringsAsFactors = FALSE)

  # Generate a scatterplot comparing the susie estimates vs. the
  # mvsusie estimates.
  pdat <- data.frame(bed     = susie$Bed_File,
                     susie   = susie$odds,
                     mvsusie = mvsusie$odds)
  scatterplots[[trait]] <-
    ggplot(pdat,aes(x = susie,y = mvsusie)) +
      geom_point() +
      geom_abline(intercept = 0,slope = 1,color = "magenta",
                  linetype = "dotted") +
      theme_cowplot(font_size = 10)
}
