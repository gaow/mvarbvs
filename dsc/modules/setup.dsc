
# Modules to provide data
# Real or simulated

# Module output
# =============
# $data: full data
# $sumstats: summary statistics

full_data: sim_utils.R + R(data =readRDS(dataset);
            data$X = as.matrix(data$X[,get_center(subset, ncol(data$X))]);
            r = round(cor(data$X), 4);
            saveRDS(r, ld_mat);
            write.table(r,ld_file,quote=F,col.names=F,row.names=F))
  tag: "full"
  dataset: Shell{head -150 ${data_file}}
  subset: NULL
  $data: data
  $top_idx: NA
  $ld_file: file(ld)
  $ld_mat: file(rds)
        
lite_data(full_data):
  tag: "2k"
  subset: 2000
             
liter_data(full_data):
  tag: "1k"
  subset: 1000
            
two_effect(full_data):
  tag: "two"
  subset: 2
             
dap_g_data(full_data): R(X = readRDS(dataset)$X;
              r = round(cor(X), 4);
              saveRDS(r, ld_mat);
              write.table(r,ld_file,quote=F,col.names=F,row.names=F)) + \
              dap_g_paper.R + R(data = list(X=X,Y=Y,true_coef=B))  
  tag: "dap_g"
  dataset: Shell{cat ${dap_g_data}}
             
get_sumstats: regression.R + R(res = mm_regression(as.matrix(data$X), 
                                                   as.matrix(data$Y), data$Z))
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
