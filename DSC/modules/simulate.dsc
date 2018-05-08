
# base_sim:
# - A base simulator of 2 independent multivariate effects
# - using MultivariateMixture
# original_Y：
# - do not simulate data, just use original

base_sim: lib_regression_simulator.py + \
                regression_simulator.py + \
                Python(data['Y'] = simulate_main(data, conf))
  data: $data
  top_idx: $top_idx
  n_signal: 3
  n_traits: 2
  eff_mode: mash_low_het
  swap_eff: raw(True)
  tag: sim1
  @ALIAS: conf = Dict(!data, !eff_mode)
  $data: data

original_Y(base_sim):
  eff_mode: original
