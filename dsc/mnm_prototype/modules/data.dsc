# Modules to provide data
# Real or simulated

# Module output
# =============
# $data: full data
# $sumstats: summary statistics


full_data: misc.R + R(data = readRDS(dataset);
            X = center_scale(filter_X(data$X[,get_center(subset, ncol(data$X))], missing_rate_cutoff, maf_cutoff)))
  tag: "full"
  dataset: Shell{head -${n_dataset} ${data_file}}
  subset: NULL
  missing_rate_cutoff: 0.05
  maf_cutoff: 0.05
  $X: X
  $Y: data$y_res
  $N: nrow(X)

lite_data(full_data):
  tag: "2k"
  subset: 2000
             
small_data(full_data):
  tag: "1k"
  subset: 1000

tiny_data(full_data):
  tag: "300"
  subset: 300