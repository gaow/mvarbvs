FLASHFMwithFINEMAP: flashfmwithFINEMAP.R
  @CONF: R_libs = (flashfm@jennasimit/flashfm)
  sumstats: $sumstats
  suffstats: $suffstats
  ld: $ld
  aaf: $aaf
  prefix: file(flashfmFINEMAP)
  $result: result
  
FLASHFMwithJAM: flashfmwithJAM.R
  @CONF: R_libs = (flashfm@jennasimit/flashfm)
  sumstats: $sumstats
  suffstats: $suffstats
  ld: $ld
  aaf: $aaf
  prefix: file(flashfmJAM)
  $result: result