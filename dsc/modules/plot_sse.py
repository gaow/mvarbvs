
def plot_sse(coef, true_coef, in_set, ld, plot_prefix):
    reg = RegressionData()
    reg.set_xcorr(ld)
    in_set = np.sum(np.array(in_set), axis = 0)
    coef = np.array(coef)
    if true_coef is not None:
        true_coef = np.array(true_coef)
    for j in range(coef.shape[1]):
        plot_file = f'{plot_prefix}.{j+1}.pdf'
        reg.plot_property_vector(coef[:,j], None,
                                 xz_cutoff = (0, 0.8), out = plot_file,
                                 conf = {'title': f'Response {j+1}', 
                                    'ylabel': 'effect size estimate', 
                                    'zlabel': 'In 95 CI set'})
