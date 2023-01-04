# R2HESS version 1.0.1
# https://www.mrc-bsu.cam.ac.uk/software/
# I obtained version 1.99 from the authors via email
mthess: misc.R + mthess.R
  @CONF: R_libs = R2HESS
  X: $X
  Y: $Y
  data_dir: file(datadir)
  $result: result