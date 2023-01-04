# Modules to provide data
# Real or simulated

# Module output
# =============
# $data: full data
# $sumstats: summary statistics


full_data: misc.R + R(data = readRDS(dataset);
            X = filter_X(data$X, missing_rate_cutoff, maf_cutoff);
            X = susieR:::set_X_attributes(X[,get_center(subset, ncol(X))]);
            X_scaled = t((t(X) - attributes(X)[["scaled:center"]]) / attributes(X)[["scaled:scale"]]);
            r = cor(X);
            var_Y = compute_cov_flash(data$y_res, flash_error))
  @CONF: R_libs = flashier@willwerscheid/flashier
  tag: "full"
  dataset: Shell{head -${n_dataset} ${data_file}}
  subset: NULL
  missing_rate_cutoff: 0.05
  maf_cutoff: 0.05
  flash_error: file(rds)
  $X: X_scaled
  $G: X
  $Y: data$y_res
  $ld: r
  $var_Y: var_Y
  $N: nrow(X)

medium_data(full_data):
  tag: "2k"
  subset: 2000

small_data(full_data):
  tag: "1k"
  subset: 1000

tiny_data(full_data):
  tag: "300"
  subset: 300
