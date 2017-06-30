import pandas as pd
import numpy as np
import os, random

def shuffle(df, axis=1):     
        return df.apply(np.random.permutation, axis=axis)

class PhenotypeSimulator:
    def __init__(self, genotype_file):
        self.gfile = genotype_file
        self.phenotype = {}
        self.beta = {}
        self.pid = None
        
    def get_genes(self, limit = 5):
        res = pd.HDFStore(self.gfile).keys()
        if len(res) > limit:
            res = res[:limit]
        return res
    
    def get_X(self, table):
        return pd.read_hdf(self.gfile, table)
    
    def permute_X(self, tables, save_to):
        if os.path.isfile(save_to):
            os.remove(save_to)
        for table in tables:
            X = pd.read_hdf(self.gfile, table)
            X = shuffle(X)
            X.to_hdf(save_to, table, mode = 'a', complevel = 9, complib = 'zlib')
        self.gfile = save_to
    
    def get_ld(self, tables, save_to = None):
        '''r^2 based LD calculation'''
        ld = {table: pd.read_hdf(self.gfile, table).transpose().corr(method = 'pearson') for table in tables}
        ld = {key: (np.power(value, 2) * np.sign(value)).astype(np.float16) for key, value in ld.items()}
        if save_to is not None:
            if os.path.isfile(save_to):
                os.remove(save_to)
            for key in ld:
                ld[key].to_hdf(save_to, key, mode = 'a', complevel = 9, complib = 'zlib')
        return ld
    
    def load_ld(self, tables, fn):
        ld = {}
        for table in tables:
            ld[table] = pd.read_hdf(fn, table)
        return ld
    
    def ld_heatmap(self, corrmat, out):
        import seaborn as sns
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        sns.heatmap(corrmat, ax = ax, vmin=-1, vmax=1, square=True, xticklabels = False, yticklabels = False)
        plt.savefig(out, dpi = 500)
        
    def generate_betamix(self, nbeta, mus, sigmas, pis, pi0 = 0):
        '''beta ~ \pi_0\delta_0 + \sum \pi_i N(0, sigma_i)
        sigma here is a nbeta list or nbeta * nbeta matrix
        '''
        if isinstance(sigmas, list):
            sigmas = np.diag(sigmas)
        assert (len(pis), len(pis)) == sigmas.shape
        masks = np.random.multinomial(1, pis, size = nbeta)
        mix = np.random.multivariate_normal(mus, sigmas, nbeta)
        return np.sum(mix * masks, axis = 1) * np.random.binomial(1, 1 - pi0, nbeta)
    
    def generate_y(self, X, beta, sigma, force = False):
        if self.pid in self.phenotype and force is not True:
            print('Name "{}" already exists. Use "force = True" to overwrite it'.format(self.pid))
            return self.phenotype[self.pid]
        assert X.shape[0] == len(beta)
        self.beta[self.pid] = beta.tolist()
        beta.reshape((len(beta),1))
        y = np.dot(X.T, beta) + np.random.normal(0, sigma, X.shape[1])
        y.reshape(len(y), 1)
        y = pd.DataFrame(data = y, columns = [self.pid], index = X.columns).transpose()
        self.phenotype[self.pid] = y
        return y
    
    def select_convoluted_snps(self, ld, cutoff1 = 0.8, cutoff2 = 10, cutoff3 = 0.01):
        '''based on LD matrix select SNPs in strong LD with other SNPs 
        yet are independent between themselves'''
        print('Count strong LD')
        strong_ld_count = ((np.absolute(ld) > cutoff1) * ld).sum(axis = 0).sort_values(ascending = False)
        strong_ld_count = strong_ld_count[strong_ld_count > cutoff2]
        print('Filter by LD')
        exclude = []
        for x in strong_ld_count.index:
            if x in exclude:
                continue
            for y in strong_ld_count.index:
                if y in exclude or y == x:
                    continue
                if np.absolute(ld[x][y]) > cutoff3:
                    exclude.append(y)
        print('Done')
        return [i for i, x in enumerate(strong_ld_count.index) if not x in exclude]
    
    def swap_beta(self, beta, strength_index):
        '''Set tops of beta to tops in strength_index'''
        nb = [0] * len(beta)
        beta = sorted(beta, key=abs, reverse=True)
        for item in strength_index:
            nb[item] = beta.pop(0)
        random.shuffle(beta)
        for idx in range(len(nb)):
            if not idx in strength_index:
                nb[idx] = beta.pop(0)
        assert len(beta) == 0
        return np.array(nb)
    
    def set_id(self, name):
        self.pid = name

        
class BetaDist:
    '''Reproducing simulated distributions of Stephens 2017 (ASH paper)'''
    def __init__(self):
        self.pi0 = 0
        self.pis = [None]
        self.mus = [None]
        self.sigmas = [None]
        
    def set_pi0(self, pi0):
        self.pi0 = pi0
        
    def set_spiky(self):
        self.pis = [0.4,0.2,0.2,0.2]
        self.mus = [0,0,0,0]
        self.sigmas = [0.25,0.5,1,2]
    
    def set_near_normal(self):
        self.pis = [2/3,1/3]
        self.mus = [0,0]
        self.sigmas = [1,2]
        
    def set_flat_top(self):
        self.pis = [1/7] * 7
        self.mus = [-1.5, -1, -.5 , 0, .5, 1, 1.5]
        self.sigmas = [0.5] * 7
        
    def set_skew(self):
        self.pis = [1/4,1/4,1/3,1/6]
        self.mus = [-2,-1,0,1]
        self.sigmas = [2,1.5,1,1]
        
    def set_big_normal(self):
        self.pis = [1]
        self.mus = [0]
        self.sigmas = [4]

    def set_bimodal(self):
        self.pis = [0.5, 0.5]
        self.mus = [-2, 2]
        self.sigmas = [1, 1]
        
    def __str__(self):
        params = ' + '.join(["{} N({}, {}^2)".format(x,y,z) for x, y, z in zip(self.pis, self.mus, self.sigmas)])
        return '{:.3f} \delta_0 + {:.3f} [{}]'.format(self.pi0, 1 - self.pi0, params)