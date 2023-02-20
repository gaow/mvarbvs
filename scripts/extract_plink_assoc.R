# TO DO: Explain here what this script does, and how to use it.
#
# Directory on gardner:
# /gpfs/data/stephens-lab/finemap-uk-biobank/data/raw/BloodCells/gwas_maf001_info6
#
library(data.table)
traits <- c("Basophill_perc","Eosinophill_perc","HLR_perc","Haemoglobin",
            "Lymphocyte_perc","MCV","MSCV","Monocyte_perc","Neutrophill_perc",
            "PDW","Platelet_count","Plateletcrit","RBC_count","RDW",
            "Reticulocyte_perc","WBC_count")
chr <- 8
pos <- c(118878922,119377965)
plink <- NULL
for (trait in traits) {
  cat(trait,"\n")
  datfile <- sprintf("bloodcells_gwas_chr%d.%s.glm.linear",chr,trait)
  dat <- fread(datfile,sep="\t",header = TRUE)
  class(dat) <- "data.frame"
  dat <- subset(dat,POS >= pos[1] & POS <= pos[2])
  plink <- cbind(plink,as.numeric(dat$P))
}
plink <- as.data.frame(plink)
rownames(plink) <- dat$ID
colnames(plink) <- traits
print(t(t(sapply(plink,min))))
