def summarize_LD(X, ld_input, ld_plot):
    data = RegressionData()
    data.X = X
    data.set_xcorr(ld_input)
    data.plot_xcorr(ld_plot)
    return data.get_representative_features()

def simulate_main(data, conf):
    '''
    data: $data
    top_idx: $top_eff
    n_signal: 3
    n_traits: 2
    eff_mode: mash_low_het
    swap_eff: raw(True)
    tag: sim1
    @ALIAS: conf = Dict(!data, !eff_mode)
    $data: data
    '''
    reg = RegressionData()
    reg.X = data['X']
    if eff_mode == 'mash_low_het':
        if n_traits < 2:
            raise ValueError(f'Cannot simulate {n_traits} under mode {eff_mode}')
        data['Y'] = mash_low_het(reg, conf)
    elif eff_mode == 'original':
        data['Y'] = original_y(data, conf)
    else:
        raise ValueError(f'Mode {eff_mode} is not implemented.')
    return data
        
def original_y(data, conf):
    return np.vstack(data['Y'].values()).T
    
def mash_low_het(reg, conf):
    return None
