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
  $plot_file: file(pdf)

plot_caviar(plot_finemap): plot_caviar.R
plot_dap(plot_finemap): plot_dap.R
