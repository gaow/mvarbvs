
# Modules to provide data
# Real or simulated

# Module output
# =============
# $data: full data
# $sumstats: summary statistics

full_data: R(data =readRDS(${data_file});
            if (end>start) data$X = as.matrix(data$X[,start:end]);
            r2 = cor(data$X);
            saveRDS(r2 ^ 2 * sign(r2), ld_mat);
            write.table(r2,ld_file,quote=F,col.names=F,row.names=F))
  tag: full
  start, end: (0, 0)
  $data: data
  $top_idx: raw(NULL)
  $ld_file: file(ld)
  $ld_mat: file(rds)
        
lite_data(full_data):
  tag: 2k
  start, end: (2500, 4500)
             
liter_data(full_data):
  tag: 1k
  start, end: (3000, 4000)           
            
two_effect(full_data):
  tag: two
  start, end: (3500, 3501)

get_sumstats: regression.R + R(res = mm_regression(as.matrix(data$X), 
                                                   as.matrix(data$Y)))
  @CONF: R_libs = abind
  data: $data
  $sumstats: res
                                                   
summarize_ld: lib_regression_simulator.py + \
                regression_simulator.py + \
                Python(res = summarize_LD(data['X'], ld_mat, ld_plot))
  data: $data
  ld_mat: $ld_mat
  $ld_plot: file(png)
  $top_idx: res
