# R2HESS version 1.0.1
# https://www.mrc-bsu.cam.ac.uk/software/
mthess: mthess.R
  @CONF: R_libs = R2HESS
  X: $X
  Y: $Y
  data_dir: file(datadir)
  $result: result

# https://github.com/hruffieux/atlasqtl
atlasqtl: atlasqtl.R
  @CONF: R_libs = R2HESS
  X: $X
  Y: $Y
  meta: $meta
  $result: result