# Modules to evaluate various methods
# for finemapping-m

# Module input
# ============
# $fit: see fit.dsc
# $posterior: see fit.dsc

# Module output
# =============
# ? an object for diagnosis

diagnose: elbo_mnm.R
  data: $data
  model: $model
  fitted: $fitted
  posterior: $posterior
  $diagnosed: elbo