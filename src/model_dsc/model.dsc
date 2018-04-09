#!/usr/bin/env dsc

get_data: R(dat = readRDS(data_file))
  # FIXME: see 20171103_MNMASH_Data.ipynb for GTEx multitissue data preparation
  # and implement it more formally here
  data_file: $(data_file)
  $data: dat

original_Y: Python(data['Y'] = numpy.vstack(data['Y'].values()).T)
  # do not simulate data, just use original
  data: $data
  $data: data

init_model: init_mnm.py
  data: $data
  (U, grid, pi): (raw({'identity':np.identity(2),'single_1':np.array([[1,0],[0,0]]),'single_2':np.array([[0,0], [0,1]]), 'all_in':np.ones((2,2))}), (0.9,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02), (0.5,1))
  V: empirical
  pi0: 0
  $data: data

fit: fit_mnm.R
  data: $data
  $model: model

diagnose: elbo_mnm.py
  model: $model
  $summary: summary

DSC:
  run:
    first_pass: get_data * original_Y * init_model * fit * diagnose
  output: mnm_model
  exec_path: modules
  global:
    data_file: ~/Documents/GTExV8/Thyroid.Lung.FMO2.filled.rds
