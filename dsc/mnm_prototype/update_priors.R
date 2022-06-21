print('artificial_mixture_ukb')
priors = readRDS('~/GitHub/mvarbvs/dsc/ukb-bloodcells/ukb_prior_simulation_20220429.rds') 
files = read.table('~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss_20220429/artificial_mixture_ukb/analysis_units.txt', stringsAsFactors = FALSE)$V1
for(i in 1:length(files)){
  print(i)
  dat = readRDS(paste0('~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss_20220429/artificial_mixture_ukb/', files[i], '.rds'))
  dat$meta$prior$ED = list(xUlist = priors[[dat$meta$eff_mode]]$ED$U, pi = priors[[dat$meta$eff_mode]]$ED$w, null_weight=0)
  # dat$meta$prior$ED_ddcan = list(xUlist = priors[[dat$meta$eff_mode]]$ED_can$U, 
  #                                pi = priors[[dat$meta$eff_mode]]$ED_can$w, null_weight=0)
  saveRDS(dat, paste0('~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss_20220429/artificial_mixture_ukb/', files[i], '.rds'))
}

print('ukb_bloodcells_mixture')
priors = readRDS('~/GitHub/mvarbvs/dsc/ukb-bloodcells/ukb_prior_simulation_20220429.rds')
files = read.table('~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss_20220429/ukb_bloodcells_mixture/analysis_units.txt', stringsAsFactors = FALSE)$V1
for(i in 1:length(files)){
  print(i)
  dat = readRDS(paste0('~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss_20220429/ukb_bloodcells_mixture/', files[i], '.rds'))
  dat$meta$prior$ED = list(xUlist = priors[[dat$meta$eff_mode]]$ED$U, pi = priors[[dat$meta$eff_mode]]$ED$w, null_weight=0)
  # dat$meta$prior$ED_ddcan = list(xUlist = priors[[dat$meta$eff_mode]]$ED_can$U, 
  #                                pi = priors[[dat$meta$eff_mode]]$ED_can$w, null_weight=0)
  saveRDS(dat, paste0('~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss_20220429/ukb_bloodcells_mixture/', files[i], '.rds'))
}

