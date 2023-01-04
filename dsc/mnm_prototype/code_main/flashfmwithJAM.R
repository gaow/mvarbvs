## FLASHFMwithJAM

library(flashfm)
library(data.table)
library(R2BGLiMS)

d1 = readRDS('output/ukb_rss_small2_20220619/data_ukb/data_ukb_1.rds')
d2 = readRDS('output/ukb_rss_small2_20220619/artificial_mixture_ukb_small2/data_ukb_1_artificial_mixture_ukb_small2_1.rds')
sumstats = d2$sumstats
suffstats = d2$suffstats
aaf = d2$aaf
ld = d1$ld
prefix='test_flashfm_JAM'

beta_input = list()
J = nrow(sumstats$bhat)
for(r in 1:ncol(sumstats$bhat)){
  beta_input[[r]] = sumstats$bhat[,r]
  names(beta_input[[r]]) = paste0('rs',1:J,'id')
}
names(beta_input) = paste0('t',1:ncol(sumstats$bhat))

# LD info can be provided either as R data objects or R RDS files
if (is.character(ld)) {
  LD = readRDS(ld)
} else {
  LD = ld
}
rownames(LD) = colnames(LD) = names(aaf) = paste0('rs',1:J,'id')

ybar = numeric(ncol(sumstats$bhat))
names(ybar) = paste0('t',1:ncol(sumstats$bhat))
covY = as.matrix(suffstats$YtY) / (suffstats$N-1)
colnames(covY) = rownames(covY) = paste0('t',1:ncol(sumstats$bhat))

result <- FLASHFMwithJAMhat(beta1=beta_input, 
                         corX = LD, 
                         raf = aaf,
                         ybar = LD,
                         N = rep(suffstats$N, ncol(sumstats$bhat)), 
                         save.path="DIRtmp",
                         TOdds=1,
                         covY = covY,
                         cpp=0.99,
                         NCORES=5) 

saveRDS(result, paste0(prefix,'.rds'))
