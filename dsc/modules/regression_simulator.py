def summarize_LD(X, ld_input, ld_plot):
    data = RegressionData()
    data.X = X
    data.set_xcorr(ld_input)
    data.plot_xcorr(ld_plot)
    return data.get_representative_features()

def simulate_main(data, c, plot_prefix):
    '''
    data: $data
    top_idx: $top_eff
    n_signal: 3
    n_traits: 2
    eff_mode: mash_low_het
    swap_eff: raw(True)
    keep_ld: raw(True)
    tag: sim1
    @ALIAS: conf = Dict(!data, !eff_mode)
    $data: data
    '''
    reg = RegressionData()
    reg.X = data['X'].astype(float)
    if eff_mode == 'original':
        c['swap_eff'] = False
    if c['swap_eff'] and c['top_idx'] is None:
        raise ValueError(f'"top_idx" variable is not set by an upstream module')
    if eff_mode == 'mash_low_het':
        if c['n_traits'] < 2:
            raise ValueError(f'Cannot simulate {c["n_traits"]} under mode {eff_mode}')
        data['true_coef'] = mash_low_het(data, reg, c)
    elif eff_mode == 'original':
        data['true_coef'] = original_y(data, reg, c)
    else:
        raise ValueError(f'Mode {eff_mode} is not implemented.')
    if c['center_data']:
        reg.center_data()
    data['X'] = reg.X
    data['Y'] = reg.Y
    if data['true_coef'] is not None:
        for j in range(data['true_coef'].shape[1]):
            plot_file = f'{plot_prefix}.{j+1}.pdf'
            reg.plot_property_vector(data['true_coef'][:,j], 
                                 [np.absolute(x)>0 for x in data['true_coef'][:,j]], 
                                 xz_cutoff = None, out = plot_file,
                                conf = {'title': f'Response {j+1}', 
                                        'ylabel': 'effect size', 'zlabel': ''})
    if data['Y'].shape[1] == 1:
        data['V'] = np.cov(data['Y'][:,0])
    else:
        data['V'] = np.cov(data['Y'], rowvar = False)
    return data
        
def original_y(data, reg, c):
    if isinstance(data['Y'], pd.DataFrame):
        reg.Y = np.vstack(data['Y'].values()).T
    else:
        reg.Y = data['Y']
    return None if 'true_coef' not in data else data['true_coef']
    
def mash_low_het(data, reg, c):
    if not c['keep_ld']:
        reg.permuate_X_columns()
        data['X'] = reg.X
    eff = MultivariateMixture((data['X'].shape[1], c['n_traits']))
    eff.set_low_het()
    eff.apply_grid()
    eff.get_effects()
    if c['swap_eff']:
        eff.swap_top_effects(c['top_idx'])
    eff.sparsify_effects(c['n_signal'])
    reg.Y = eff.get_y(reg, ResidualVariance(c['residual_mode']).apply(eff))
    return eff.coef
