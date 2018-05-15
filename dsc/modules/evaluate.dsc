# Modules to evaluate various methods
# for finemapping-m

# Module input
# ============
# $fit: see fit.dsc
# $result: see fit.dsc

# Module output
# =============
# ? an object or plot for diagnosis

plot_finemap: plot_finemap.R
  @CONF: R_libs = (dplyr, ggplot2, cowplot)
  result: $posterior
  top_rank: 10
  $plot_file: file(png)

plot_caviar(plot_finemap): plot_caviar.R
plot_dap(plot_finemap): plot_dap.R

plot_sse: lib_regression_simulator.py + \
            plot_sse.py + \
            Python(plot_sse(result['PosteriorMean'], data['true_coef'],
                            result['in_CI'], ld_mat, plot_file))
  @CONF: python_modules = seaborn
  data: $data
  result: $posterior
  ld_mat: $ld_mat
  $plot_file: file(plot_file)

plot_susie: plot_susie.py + Python(purity = plot_sets(result['in_CI'], 
                                             result['n_in_CI'],
                                             result['set_lfsr'],
                                             data['true_coef'],
                                             ld_mat,
                                             seg))
  @CONF: python_modules = seaborn
  data: $data
  result: $posterior
  ld_mat: $ld_mat
  $seg: file(seg)
  $purity: purity
