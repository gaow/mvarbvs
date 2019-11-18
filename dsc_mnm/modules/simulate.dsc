# A module to prepare for simulation configurations
oracle_generator: lib_regression_simulator.py + \
                regression_simulator.py + \
                Python(configurations = get_config(effects, n_traits, ${grid}, dict(mixture_1=${mixture_1})))
  effects: ('identity', 'low_het', 'mid_het', 'high_het', 'shared', 'singleton', 'singleton_1', 'mixture_1')
  $configurations: configurations
  n_traits: (${R})

# MASH 'identity' matrix
# No correlated effects across conditions
identity: lib_regression_simulator.py + \
                regression_simulator.py + \
                Python(res = simulate_main(dict(X=X,Y=Y), conf))
  @ALIAS: conf = Dict(!X, !Y)
  @CONF: python_modules = seaborn
  X: $X
  Y: $Y
  n_traits: ${R}
  # set signal to <0 to use a default setting
  n_signal: ${C}
  eff_mode: "identity"
  eff_grid: ${grid}
  residual_mode: "identity"
  keep_ld: True
  center_data: True
  # per-condition PVE when evaluated in univeriate framework
  pve: ${pve}
  $Y: res['Y']
  $R: res['Y'].shape[1]
  $J: res['X'].shape[1]
  $pve_out: conf['pve']
  $meta: dict(true_coef=res['true_coef'], residual_variance=res['residual_variance'], original_Y=Y)

# MASH 'simple_het_1' matrix
# where off-diagonal is 0.8 sharing
low_het(identity):
  eff_mode: "low_het" 

# MASH 'simple_het_2' matrix
# where off-diagonal is 0.5 sharing
mid_het(identity):
  eff_mode: "mid_het" 

# MASH 'simple_het_3' matrix
# where off-diagonal is 0.25 sharing
high_het(identity):
  eff_mode: "high_het" 

# MASH 'equal_effects' matrix
shared(identity):
  eff_mode: "shared"

# MASH 'singleton' matrix
# where we randomly let one of the conditions to have non-zero equal_effect
# and other conditions have zero effects
singleton(identity):
  eff_mode: "singleton"

# MASH 'singleton' matrix
# where we fix one of the conditions (in this case, the first condition) to have non-zero equal_effect
# and other conditions have zero effects
singleton_first(identity):
  eff_mode: "singleton_1"

# MASH manual mixture setup
mixture01(identity):
  eff_mode: "mixture_1"
  # it is okay if these weights do not sum to one because they will be normalized down the road
  mixture_dist: ${mixture_1}