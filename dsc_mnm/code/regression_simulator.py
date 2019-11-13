import numpy as np

def summarize_LD(X, ld_input, ld_plot):
    data = RegressionData()
    data.X = X
    data.set_xcorr(ld_input)
    data.plot_xcorr(ld_plot)
    return data.get_representative_features()

def get_n_signal():
    n = [1,2,3,4,5,6,7,8,9,10]
    pi = [0.5,0.3,0.2,0,0,0,0,0,0,0]
    one_hot = list(np.random.multinomial(1, pi))
    return n[one_hot.index(1)]

def simulate_main(data, c):
    '''
    data: $data
    n_traits: 2
    eff_mode: low_het
    keep_ld: True
    @ALIAS: conf = Dict(!data, !eff_mode)
    $data: data
    '''
    eff_mode = c['eff_mode']
    if not c['keep_ld']:
        data['X'].permuate_X_columns()
    data['n_signal'] = get_n_signal() if (not 'n_signal' in c or c['n_signal'] <=0) else int(c['n_signal'])
    reg = RegressionData()
    reg.X = data['X'].astype(float)
    if eff_mode == 'original':
        data['true_coef'] = original_y(data, reg, c)
        data['residual_variance'] = None
    elif eff_mode == 'simple_lm':
        data['true_coef'], data['residual_variance'], data['prior'] = simple_lm(data, reg, c)
    else:
        if c['n_traits'] < 2:
            raise ValueError(f'Cannot simulate {c["n_traits"]} under mode {eff_mode}')
        data['true_coef'], data['residual_variance'], data['residual_correlation'], data['prior'] = mash_sim(data, reg, c, eff_mode)
    if c['center_data']:
        reg.center_data()
    data['X'] = reg.X
    data['Y'] = reg.Y
    if data['Y'].shape[1] == 1:
        data['varY'] = np.cov(data['Y'][:,0])
    else:
        data['varY'] = np.cov(data['Y'], rowvar = False)
    return data
        
def original_y(data, reg, c):
    '''
    Use the original Y input instead of simulating any
    '''
    if 'y' in data:
        reg.Y = data.pop('y')
    elif isinstance(data['Y'], pd.DataFrame):
        reg.Y = np.vstack(data['Y'].values()).T
    else:
        reg.Y = data['Y']
    return None if 'true_coef' not in data else np.array(data['true_coef'])

def simple_lm(data, reg, c):
    '''
    simple linear model simulation
    '''
    if 'Z' in data:
        del data['Z']
    if 'y' in data:
        del data['y']
    eff = UnivariateMixture(reg.X.shape[1])
    eff.set_vanilla(c['amplitude'])
    Y = []
    coef = []
    sigma = []
    for i in range(c['n_traits']):
        eff.get_effects()
        eff.sparsify_effects(data['n_signal'])
        Y.append(eff.get_y(reg, pve = c['pve']))
        coef.append(eff.coef)
        sigma.append(eff.residual_variance)
    reg.Y = np.asmatrix(Y).T
    return np.array(coef).T, np.array(sigma), c['pve']
    
def mash_sim(data, reg, c, mode):
    if 'y' in data:
        del data['y']
    #eff = MultivariateMixtureEqualEff((data['X'].shape[1], c['n_traits']))
    eff = MultivariateMixture((data['X'].shape[1], c['n_traits']))
    eff.set_grid(c['eff_grid'])
    if mode == 'low_het':
        eff.set_low_het()
    elif mode == 'mid_het':
        eff.set_mid_het()
    elif mode == 'high_het':
        eff.set_high_het()
    elif mode == 'shared':
        eff.set_shared()
    elif mode == 'identity':
        eff.set_indep()
    elif mode.startswith('singleton'):
        eff.set_singleton(mode)
    elif mode.startswith("mixture_"):
        eff.set_manual_mixture(c['mixture_dist'])
    elif mode.startswith("group_"):
        eff.set_manual_blocks(c['n_groups'], c['group_het'], c['group_prop'])
    else:
        raise ValueError(f"Unknown mash simulation mode `{mode}`")
    eff.apply_grid()
    eff.get_effects()
    eff.sparsify_effects(data['n_signal'])
    residual_correlation = ResidualCorrelation(c['residual_mode'], eff.R).apply()
    reg.Y = eff.get_y(reg, c['pve']/data['n_signal'], residual_correlation, is_pve_per_variable = True)
    return eff.coef, eff.residual_variance, residual_correlation, eff.get_prior()

def get_config(effects, R, eff_grid, mixtures):
    res = dict()
    for mode in effects:
        eff = MultivariateMixtureEqualEff((0, R))
        eff.set_grid(eff_grid)
        if mode == 'low_het':
            eff.set_low_het()
        elif mode == 'mid_het':
            eff.set_mid_het()
        elif mode == 'high_het':
            eff.set_high_het()
        elif mode == 'shared':
            eff.set_shared()
        elif mode == 'identity':
            eff.set_indep()
        elif mode.startswith('singleton'):
            eff.set_singleton(mode)
        elif mode.startswith("mixture_"):
            eff.set_manual_mixture(mixtures[mode])
        elif mode.startswith("group_"):
            eff.set_manual_blocks(mixtures[mode]['n_groups'], mixtures[mode]['group_het'], mixtures[mode]['group_prop'])
        else:
            raise ValueError(f"Unknown mash simulation mode `{mode}`")
        eff.apply_grid()
        res[mode] = eff.get_prior()
    return res