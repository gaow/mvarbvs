# TO DO: Explain here what this script is for, and how to use it.

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

# Select SNPs with a PIP of 0.9 or more.
susie   <- subset(susie,PIP > 0.9)
mvsusie <- subset(mvsusie,PIP > 0.9)

# Add a "susie" column to the "mvsusie" data frame,
mvsusie$susie <- "no"
rows <- which(is.element(mvsusie$ID,susie$ID))
mvsusie[rows,"susie"] <- "yes"
mvsusie$susie <- factor(mvsusie$susie)

# Reorganize the columns and rows in the "mvsusie" data frame a bit.
mvsusie <- mvsusie[,c("CHR","POS","ID","REF","ALT","maf","PIP","susie",
                      "Region","CS_trait")]
names(mvsusie) <- c("chr","pos","id","ref","alt","maf","pip","susie",
                    "region","traits")
rows <- order(mvsusie$chr,mvsusie$pos)
mvsusie <- mvsusie[rows,]
rownames(mvsusie) <- NULL

# Write the "mvsusie" data frame to a CSV file.
mvsusie <- transform(mvsusie,
                     maf = round(maf,digits = 4),
                     pip = round(pip,digits = 4))
write.csv(mvsusie,"blood_cell_traits_mvsusie_sentinels.csv",
          quote = FALSE,row.names = FALSE)
