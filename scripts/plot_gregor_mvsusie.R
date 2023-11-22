# TO DO: Explain here what this script is for, and how to use it.
library(ggplot2)
library(cowplot)
library(ggrepel)
traits <- c("RBC","HGB","MCV","RDW","MSCV",          # mature red blood cell
            "Reticulocyte_perc","HLR_perc",          # immature red blood cell
            "Platelet_count","Plateletcrit","PDW",   # platelet
            "WBC","Lymphocyte_perc","Monocyte_perc", # white blood cell
            "Neutrophill_perc","Eosinophill_perc","Basophill_perc",
            "crosstrait")

# Repeat for each trait.
n <- length(traits)
scatterplots <- vector("list",n)
names(scatterplots) <- traits
mvsusie_odds_all <- NULL
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
  susie   <- subset(susie,is.finite(odds))
  mvsusie <- subset(mvsusie,is.finite(odds))

  # Include the mvsusie results in the "mvsusie_odds_all" data frame.
  beds <- mvsusie$Bed_File
  mvsusie_odds_all <- rbind(mvsusie_odds_all,
                            data.frame(bed   = mvsusie$Bed_File,
                                       trait = trait,
                                       odds  = mvsusie$odds,
                                       stringsAsFactors = FALSE))
  
  # Generate a scatterplot comparing the susie estimates vs. the
  # mvsusie estimates.
  pdat <- data.frame(bed     = susie$Bed_File,
                     susie   = susie$odds,
                     mvsusie = mvsusie$odds,
                     stringsAsFactors = FALSE)
  x <- c(pdat$susie,pdat$mvsusie)
  rows <- which(with(pdat,!(mvsusie > 5 & mvsusie/susie > 2)))
  pdat[rows,"bed"] <- ""
  scatterplots[[trait]] <-
    ggplot(pdat,aes(x = susie,y = mvsusie,label = bed)) +
      geom_point(shape = 20,color = "black",size = 1) +
      geom_abline(intercept = 0,slope = 1,color = "magenta",
                  linetype = "dotted") +
      geom_text_repel(color = "royalblue",max.overlaps = Inf,size = 1.75) +
      xlim(c(0,max(x))) +
      ylim(c(0,max(x))) +
      ggtitle(trait) +
      theme_cowplot(font_size = 8) +
          theme(plot.title = element_text(face = "plain",size = 8))
}

# Show the scatterplots and generate a PDF of the scatterplots.
# print(do.call("plot_grid",c(scatterplots,list(nrow = 5,ncol = 4))))
# TO DO.

# Here is a plot to compare the mvSuSiE enrichments across traits
# *and* annotations.
mvsusie_odds_all <- transform(mvsusie_odds_all,
                              bed   = factor(bed,beds),
                              trait = factor(trait,traits))
p <- ggplot(mvsusie_odds_all,aes(x = trait,y = bed,color = odds)) +
  geom_point(shape = 20) +
  theme_cowplot(font_size = 9)
print(p)

res <- rep(0,16)
names(res) <- traits[-17]
for (i in traits[-17]) {
  x <- subset(mvsusie_odds_all,trait == i)
  y <- subset(mvsusie_odds_all,trait == "crosstrait")
  res[i] <- cor(x$odds,y$odds)
}

hist(x$odds,n = 64,xlab = "trait/crosstrait")
