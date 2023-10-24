# Script to generate text files containing information about the
# mvSuSiE 1-SNP CSs.
library(readxl)
library(ggplot2)
library(cowplot)

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
vuckovic <- read_excel("../data/vuckovic_table_s5.xlsx",col_names = TRUE)
class(vuckovic) <- "data.frame"

# Read in the Ulirsch et al fine-mapping results (Supplementary Table
# 1 of Ulirsch et al 2019).
ulirsch <- suppressWarnings(
  read_excel("../data/ulirsch_supplementary_table_1.xlsx",
             col_names = TRUE,skip = 2))
class(ulirsch) <- "data.frame"

# Select SNPs with a PIP of 0.95 or more.
mvsusie <- subset(mvsusie,PIP >= 0.95)

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

# Add an "ulirsch_pip" column to the "mvsusie" data frame.
ulirsch <- subset(ulirsch,rsID != "NA")
ulirsch <- transform(ulirsch,rsID = factor(rsID))
pip <- tapply(ulirsch$PP,ulirsch$rsID,max)
ids <- names(pip)
mvsusie$ulirsch_pip <- as.numeric(NA)
rows1 <- which(is.element(ids,mvsusie$ID))
rows2 <- match(ids[rows1],mvsusie$ID)
mvsusie[rows2,"ulirsch_pip"] <- pip[rows1]

# Reorganize the columns and rows in the "mvsusie" data frame a bit.
mvsusie <- mvsusie[,c("CHR","POS","ID","REF","ALT","maf","PIP","susie_pip",
                      "vuckovic_pip","ulirsch_pip","Region","CS_trait")]
names(mvsusie) <- c("chr","pos","id","ref","alt","maf","pip","susie_pip",
                    "vuckovic_pip","ulirsch_pip","region","traits")
rows <- order(mvsusie$chr,mvsusie$pos)
mvsusie <- mvsusie[rows,]
rownames(mvsusie) <- NULL

# Summarize the overlap of mvsusie 1-SNP CSs vs. the other
# fine-mapping results.
susie_cs1snp <- subset(susie,PIP >= 0.95)[,"ID"]
susie_cs1snp <- unique(susie_cs1snp)
cat("susie PIPs > 0.95:",length(susie_cs1snp),"\n")
cat("mvsusie PIPs > 0.95:",nrow(mvsusie),"\n")
cat("Vuckovic et al PIPs > 0.95 (overlapping with mvsusie PIPs > 0.95):",
    sum(mvsusie$vuckovic_pip >= 0.95,na.rm = TRUE),"\n")
cat("Ulirsch et al PIPs > 0.95 (overlapping with mvsusie PIPs > 0.95):",
    sum(mvsusie$ulirsch_pip >= 0.95,na.rm = TRUE),"\n")

# Generate text files containing the 1-CS SNPs.
write.table(data.frame(susie_cs1snp),"susie_cs1snp.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(mvsusie["id"],"mvsusie_cs1snp.txt",quote = FALSE,
            row.names = FALSE,col.names = FALSE)

# Write the "mvsusie" data frame to a CSV file.
mvsusie <- transform(mvsusie,
                     maf          = format(maf,digits = 4),
                     pip          = round(pip,digits = 4),
                     susie_pip    = round(susie_pip,digits = 4),
                     vuckovic_pip = round(vuckovic_pip,digits = 4),
                     ulirsch_pip  = round(ulirsch_pip,digits = 4))
write.csv(mvsusie,"blood_cell_traits_mvsusie_cs1snp.csv",
          quote = FALSE,row.names = FALSE)
