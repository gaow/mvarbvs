# Modules to provide data
# get ukb data
data_ukb: data_ukb.R
  tag: 'full'
  dataset: Shell{head -${n_dataset} ${data_file}}
  genotype_dir: ${genotype_dir}
  ld_dir: ${ld_dir}
  varY_file: ${varY_file}
  $X: geno.file
  $Y: 0
  $ld: ld.file
  $var_Y: varY
  
