priors = readRDS('~/GitHub/mvarbvs/dsc/ukb-bloodcells/ukb_prior_simulation.rds') 
files = read.table('~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss/artificial_mixture_ukb_sim.txt', stringsAsFactors = FALSE)$V1
print('artificial_mixture_ukb')
for(i in 1:length(files)){
  print(i)
  dat = readRDS(paste0('~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss/artificial_mixture_ukb/', files[i]))
  dat$meta$prior$ED = list(xUlist = priors[[dat$meta$eff_mode]]$ED$U, pi = priors[[dat$meta$eff_mode]]$ED$w, null_weight=0)
  saveRDS(dat, paste0('~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss/artificial_mixture_ukb/', files[i]))
}

print('ukb_bloodcells_mixture')
priors = readRDS('~/GitHub/mvarbvs/dsc/ukb-bloodcells/ukb_prior_simulation.rds')
files = read.table('~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss/ukb_bloodcells_mixture_sim.txt', stringsAsFactors = FALSE)$V1
for(i in 1:length(files)){
  print(i)
  dat = readRDS(paste0('~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss/ukb_bloodcells_mixture/', files[i]))
  dat$meta$prior$ED = list(xUlist = priors[[dat$meta$eff_mode]]$ED$U, pi = priors[[dat$meta$eff_mode]]$ED$w, null_weight=0)
  saveRDS(dat, paste0('~/GitHub/mvarbvs/dsc/mnm_prototype/output/ukb_rss/ukb_bloodcells_mixture/', files[i]))
}

