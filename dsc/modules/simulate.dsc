
# base_sim:
# - A base simulator of 2 independent multivariate effects
# - using MultivariateMixture
# original_Y：
# - do not simulate data, just use original

base_sim: lib_regression_simulator.py + \
                regression_simulator.py + \
                Python(data = simulate_main(data, conf, conf['cache']))
  @CONF: python_modules = (seaborn, matplotlib, pprint)
  data: $data
  top_idx: $top_idx
  n_signal: 3
  n_traits: 2
  eff_mode: mash_low_het
  residual_mode: identity
  swap_eff: raw(True)
  keep_ld: raw(True)
  center_data: raw(True)
  cache: file(sim)
  tag: sim1
  @ALIAS: conf = Dict(!data, !eff_mode)
  $data: data
  $V: data['V']
  $N: data['Y'].shape[0]

simple_lm(base_sim):
  eff_mode: simple_lm
  amplitude: 0.5

original_Y(base_sim):
  eff_mode: original
