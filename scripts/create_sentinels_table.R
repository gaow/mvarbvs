# TO DO: Explain here what this script is for, and how to use it.
library(readxl)

# Load the susie_rss and mvsusie_rss blood cell trait fine-mapping results.
susie_outfile   <- file.path("../output/blood_cell_traits",
                             "susierss.notrem.CS_purity0.5.summary.csv.gz")
mvsusie_outfile <- paste0("../output/blood_cell_traits/",
                          "LDoriginal.Ycor.mvsusierss.CS_purity0.5.",
                          "CS_lfsr0.01.summary.csv.gz")
susie   <- read.csv(susie_outfile,header = TRUE,stringsAsFactors = FALSE)
mvsusie <- read.csv(mvsusie_outfile,header = TRUE,stringsAsFactors = FALSE)
mvsusie <- transform(mvsusie,
                     CHR    = factor(CHR,1:22),
                     Region = factor(Region),
                     REF    = factor(REF),
                     ALT    = factor(ALT))

# Read in the Vuckovic et al fine-mapping results (Supplementary Table
# S5 of Vuckovic et al 2020).
vuckovic <- read_excel("../output/blood_cell_traits/vuckovic_table_s5.xlsx")
class(vuckovic) <- "data.frame"

# Select SNPs with a PIP of 0.9 or more.
mvsusie <- subset(mvsusie,PIP > 0.9)

# Add a "susie_pip" column to the "mvsusie" data frame,
susie <- transform(susie,ID = factor(ID))
pip <- tapply(susie$PIP,susie$ID,max)
ids <- names(pip)
mvsusie$susie_pip <- as.numeric(NA)
rows1 <- which(is.element(ids,mvsusie$ID))
rows2 <- match(ids[rows1],mvsusie$ID)
mvsusie[rows2,"susie_pip"] <- pip[rows1]

# Add a "vuckovic_pip" column to the "mvsusie" data frame.
vuckovic <- transform(vuckovic,rsID = factor(rsID))
pip <- tapply(vuckovic$FM_PP,vuckovic$rsID,max)
ids <- names(pip)
mvsusie$vuckovic_pip <- as.numeric(NA)
rows1 <- which(is.element(ids,mvsusie$ID))
rows2 <- match(ids[rows1],mvsusie$ID)
mvsusie[rows2,"vuckovic_pip"] <- pip[rows1]

# Reorganize the columns and rows in the "mvsusie" data frame a bit.
mvsusie <- mvsusie[,c("CHR","POS","ID","REF","ALT","maf","PIP","susie_pip",
                      "vuckovic_pip","Region","CS_trait")]
names(mvsusie) <- c("chr","pos","id","ref","alt","maf","pip","susie_pip",
                    "vuckovic_pip","region","traits")
rows <- order(mvsusie$chr,mvsusie$pos)
mvsusie <- mvsusie[rows,]
rownames(mvsusie) <- NULL
    
# Write the "mvsusie" data frame to a CSV file.
mvsusie <- transform(mvsusie,
                     maf          = format(maf,digits = 4),
                     pip          = round(pip,digits = 4),
                     susie_pip    = round(susie_pip,digits = 4),
                     vuckovic_pip = round(vuckovic_pip,digits = 4))
write.csv(mvsusie,"blood_cell_traits_mvsusie_sentinels.csv",
          quote = FALSE,row.names = FALSE)
