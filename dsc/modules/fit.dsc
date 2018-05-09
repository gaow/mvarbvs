# workhorse(s) for finemapping

# Module input
# ============
# $data: full data; or
# $sumstats: summary statistics; or / and
# $ld: LD information

# Module output
# =============
# $fitted: for diagnostics
# $posterior: for inference

init_mnm: init_mnm.R
  # mashr comes from `dev` branch on github
  @CONF: R_libs = mashr
  V: $V
  reg: $sumstats
  # FIXME: these quantities are to be computed seperately and globally using mashr procedure
  # See http://stephenslab.github.io/gtex-eqtls/analysis/20171002_MASH_V8.html
  Sigma: empirical
  (U, grid, p): (auto, (0.9,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02), auto)
  $model: model
  $V: V

fit_mnm_debug: regression.R + elbo_mnm.R + fit_mnm.R
  @CONF: R_libs = mashr
  maxL: 5
  maxI: 20
  get_elbo: raw(T)
  data: $data
  model: $model
  V: $V
  $fitted: fitted_track
  $posterior: posterior

fit_mnm(fit_mnm_debug):
  maxI: 10
  get_elbo: raw(F)

fit_susie: fit_susie.R
  # Prior variance of nonzero effects.
  @CONF: R_libs = susieR@stephenslab/susieR
  maxL: 5
  maxI: 50
  data: $data
  $posterior: posterior
  $fitted: fitted

fit_varbvs(fit_susie): setup_varbvs.R + fit_varbvs.R
  @CONF: R_libs = varbvs@pcarbo/varbvs/varbvs-R
  sa: 1

fit_caviar: fit_caviar.R + \
             R(posterior = finemap_mcaviar(sumstats[1,,]/sumstats[2,,], 
                                            ld, args, prefix=cache))
  @CONF: R_libs = (dplyr, magrittr)
  sumstats: $sumstats
  ld: $ld_file
  args: -c 1, -c 2
  cache: file(CAVIAR)
  $posterior: posterior

fit_finemap(fit_caviar): fit_finemap.R + \
             R(posterior = finemap_mvar(sumstats[1,,]/sumstats[2,,], 
                                        ld, N, k,
                                        args, prefix=cache))
  N: $N
  k: R(rep(1/5,5)), (0.6,0.25,0.1,0.05)
  args: --regions 1 --prior-std 0.4 --n-causal-max 5
  cache: file(FM)

fit_dap: fit_dap.py + Python(posterior = dap_batch(data['X'], data['Y'], cache, args))
  data: $data
  args: -ld_control 0.25
  cache: file(DAP)
  $posterior: posterior

fit_dap_z: fit_dap.py + Python(posterior = dap_batch_z(sumstats[0,:,:]/sumstats[1,:,:],
                                                       ld, cache, args))
  sumstats: $sumstats
  ld: $ld_file
  args: -t 4
  cache: file(DAP)
  $posterior: posterior

# fit_dap_mv(fit_dap): fit_dap.py + Python(res = dap_mv())
# fit_dap_mv_ss(fit_dap): fit_dap.py + Python(res = dap_mv_ss())
