def summarize_LD(X, ld_input, ld_plot):
    data = RegressionData()
    data.X = X
    data.set_xcorr(ld_input)
    data.plot_xcorr(ld_plot)
    return data.get_representative_features()
