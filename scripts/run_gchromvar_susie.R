# Perform gchromVAR analysis of SuSiE fine-mapping results for UK
# Biobank blood cell traits.
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

# Load the susie_rss results. 

susie_outfile   <- file.path("../output/blood_cell_traits",
                             "susierss.notrem.CS_purity0.5.summary.csv.gz")
susie <- read.csv(susie_outfile,header = TRUE,stringsAsFactors = FALSE)
susie <- transform(susie,
                   CHR    = factor(CHR,1:22),
                   Region = factor(Region),
                   trait  = factor(trait),
                   REF    = factor(REF),
                   ALT    = factor(ALT))

# Limit to results with PIP > 0.01.
susie <- subset(susie,PIP > 0.01)

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
traits    <- levels(susie$trait)
bed_files <- paste(traits,"bed",sep = ".")
for (i in traits) {
  dat  <- subset(susie,trait == i)
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
susie_se   <- importBedScore(rowRanges(hema_dat),bed_files,colidx = 5)
wdev_susie <- computeWeightedDeviations(hema_dat,susie_se)

# Extract the z-scores and compute p-values.
gchromvar_susie <- assays(wdev_susie)$z
names(dimnames(gchromvar_susie)) <- c("trait","celltype")
gchromvar_susie <- melt(gchromvar_susie)
names(gchromvar_susie) <- c("trait","celltype","zscore")
gchromvar_susie$pval <-
  -log10(pnorm(gchromvar_susie$zscore,lower.tail = FALSE))

# Save the results as a CSV file.
dat <- transform(gchromvar_susie,
                 zscore = round(zscore,digits = 4),
                 pval   = round(pval,digits = 4))
write.csv(dat,"gchromvar_susie.csv",quote = FALSE,row.names = FALSE)
