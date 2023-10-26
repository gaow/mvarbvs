# Perform a gchromVAR analysis of the mvSuSiE fine-mapping results for
# UK Biobank blood cell traits.
library(chromVAR)
library(gchromVAR)
library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)

# For reproducibility, set the seed.
set.seed(1)

# Load the mvsusie_rss results. 
mvsusie_outfile <- paste0("../output/blood_cell_traits/",
                          "LDoriginal.Ycor.mvsusierss.CS_purity0.5.",
                          "CS_lfsr0.01.summary.csv.gz")
mvsusie <- read.csv(mvsusie_outfile,header = TRUE,stringsAsFactors = FALSE)
mvsusie <- transform(mvsusie,
                     CHR    = factor(CHR,1:22),
                     Region = factor(Region),
                     REF    = factor(REF),
                     ALT    = factor(ALT))

# Limit to results with PIP > 0.01.
mvsusie <- subset(mvsusie,PIP > 0.01)

# Build a set of consensus, equal-width peaks spanning the hematopoesis
# phenotypes.
counts_file <- paste0(system.file("extdata/paper",package = "gchromVAR"),
                      "/hemeCounts.tsv.gz")
peaks_file <- paste0(system.file("extdata/paper",package = "gchromVAR"),
                     "/hemePeaks.bed.gz")
peaksdf <- data.frame(fread(paste0("zcat < ",peaks_file)))
peaks <- makeGRangesFromDataFrame(peaksdf,
                                  seqnames = "V1",
                                  start.field = "V2",
                                  end.field = "V3")
counts <- data.matrix(fread(paste0("zcat < ", counts_file)))
hema_dat <- SummarizedExperiment(assays = list(counts = counts),
                                 rowData = peaks, 
                                 colData = DataFrame(names = colnames(counts)))
hema_dat <- addGCBias(hema_dat,genome = BSgenome.Hsapiens.UCSC.hg19)


## Encode the fine-mapping results for each trait as a BED file.
traits <- c("Basophill_perc","Eosinophill_perc","Haemoglobin","HLR_perc",
            "Lymphocyte_perc","MCV","Monocyte_perc","MSCV","Neutrophill_perc",
            "PDW","Platelet_count","Plateletcrit","RBC_count","RDW",
            "Reticulocyte_perc","WBC_count")
cs_traits <- lapply(strsplit(mvsusie$CS_trait,"|",fixed = TRUE),trimws)
bed_files <- paste(traits,"bed",sep = ".")
for (i in traits) {
  rows <- which(sapply(cs_traits,function (x) any(is.element(i,x))))
  dat  <- mvsusie[rows,]
  dat  <- dat[c("CHR","POS","Region","PIP")]
  names(dat) <- c("chr","start","region","pip")
  rows <- order(dat$chr,dat$start)
  dat  <- dat[rows,]
  dat  <- transform(dat,
                    chr = paste0("chr",chr),
                    end = start + 1)
  dat     <- dat[c("chr","start","end","region","pip")]
  dat$pip <- format(dat$pip,digits = 4) 
  write.table(dat,paste(i,"bed",sep = "."),quote = FALSE,sep = "\t",
              row.names = FALSE,col.names = FALSE)
}

# Now we are ready to compute the weighted deviations using gchromVAR.
mvsusie_se   <- importBedScore(rowRanges(hema_dat),bed_files,colidx = 5)
wdev_mvsusie <- computeWeightedDeviations(hema_dat,mvsusie_se)

# Extract the z-scores and compute p-values.
gchromvar_mvsusie <- assays(wdev_mvsusie)$z
names(dimnames(gchromvar_mvsusie)) <- c("trait","celltype")
gchromvar_mvsusie <- melt(gchromvar_mvsusie)
names(gchromvar_mvsusie) <- c("trait","celltype","zscore")
gchromvar_mvsusie$pval <-
  -log10(pnorm(gchromvar_mvsusie$zscore,lower.tail = FALSE))

# Save the results as a CSV file.
dat <- transform(gchromvar_mvsusie,
                 zscore = round(zscore,digits = 4),
                 pval   = round(pval,digits = 4))
write.csv(dat,"gchromvar_mvsusie.csv",quote = FALSE,row.names = FALSE)
