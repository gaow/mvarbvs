## BayesSUR
source('code/misc.R')
d1 = readRDS('output/ukb_rss_small3_20220619/data_ukb/data_ukb_1.rds')
d2 = readRDS('output/ukb_rss_small3_20220619/artificial_mixture_ukb_small3/data_ukb_1_artificial_mixture_ukb_small3_1.rds')

X = get_genotype(d1$X)

library(BayesSUR)
m = BayesSUR(Y = as.matrix(d2$Y), X = as.matrix(X), 
             outFilePath = 'BayesSURtest',
             gammaPrior = 'hotspot',
             covariancePrior = 'HIW')
saveRDS(m, 'BayesSURtest_3.rds')
