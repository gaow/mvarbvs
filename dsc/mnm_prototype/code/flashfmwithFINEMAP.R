## FLASHFMwithFINEMAP

library(flashfm)
library(data.table)

gwas.list = list()
J = nrow(sumstats$bhat)
for(r in 1:ncol(sumstats$bhat)){
  gwas.list[[r]] = data.frame(rsid = paste0('rs',1:J,'id'),
                              chromosome = 1,
                              position = 1:J,
                              allele1 = 'A',
                              allele2 = 'C',
                              maf = pmin(aaf, 1-aaf),
                              beta = sumstats$bhat[,r],
                              se = sumstats$sbhat[,r])
}
names(gwas.list) = paste0('t',1:ncol(sumstats$bhat))

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

result = FLASHFMwithFINEMAP(gwas.list = gwas.list,
                            corX = LD,
                            raf = aaf,
                            ybar = ybar,
                            N = rep(suffstats$N, ncol(sumstats$bhat)),
                            fstub = prefix,
                            TOdds=1,
                            covY = covY,
                            cpp=.99,
                            NCORES=4,
                            FMpath = 'code/finemap_v1.4.1')
