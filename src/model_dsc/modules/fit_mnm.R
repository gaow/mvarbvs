source("utils.R")
reg = univariate_regression(data$X, data$Y)
mash_data = mashr::mash_set_data(as.matrix(reg$betahat), Shat = as.matrix(reg$sebetahat), V = as.matrix(data$V))
posterior = mashr::mash_compute_posterior_matrices(data$fitted_g, mash_data)
model = posterior
                                        # model = {'post_mean_mat' : np.zeros(R),
                                        #         'post_mean2_mat' : None
                                        #         'neg_prob_mat' : None
                                        #         'zero_prob_mat' : None
                                        #         'is_common_cov' : None
                                        #         'Sigma : Sigma
                                        #         'V : None
                                        #         'pi : None
                                        #         'posterior_weights : None
                                        #         'grid : None
                                        #         'l10bf : None
                                        #         'lik : {'relative_likelihood' : None,
                                        #                          'lfactor': None,
                                        #                          'marginal_loglik': None,
                                        #                          'loglik': None,
                                        #                          'null_loglik': None,
                                        #                          'alt_loglik': None}
