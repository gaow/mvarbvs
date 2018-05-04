# Modules to provide data
# Real or simulated

# Module output
# =============
# $data: full data
# $sumstats: summary statistics

get_data: Shell(ln -sf `realpath ${data_file}` $data)
  # FIXME: see 20171103_MNMASH_Data.ipynb for GTEx multitissue data preparation
  # and implement it more formally here
  $data: file(rds)

original_Y: Python(data['Y'] = numpy.vstack(data['Y'].values()).T)
  # do not simulate data, just use original
  data: $data
  $data: data

get_sumstats: regression.R + R(res = mm_regression(data$X, data$Y); r2 = cor(data$X)^2)
  @CONF: R_libs = abind
  data: $data
  $sumstats: res
  $ld: r2