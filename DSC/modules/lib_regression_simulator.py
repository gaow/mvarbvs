import numpy as np
import os, copy
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pprint import pformat
from dsc.utils import load_rds

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError:
            raise AttributeError(item)

    def __deepcopy__(self, memo):
        return dotdict(copy.deepcopy(dict(self)))
    
class RegressionData(dotdict):
    def __init__(self, X = None, Y = None, Z = None):
        # FIXME: check if inputs are indeed numpy arrays
        self.debug = dotdict()
        self.x_centered = self.y_centered = self.z_centered = False
        self.X = None
        self.Y = None
        self.Z = None
        self.xcorr = None

    def get_summary_stats(self):
        '''
        Computer univariate regression for every X_j (N by 1) and Y_r (N by 1)
        Bhat: J by R matrix of estimated effects
        Shat: J by R matrix of SE of Bhat
        '''
        if self.Z is not None:
            self.remove_covariates()
        # Compute betahat
        XtX_vec = np.einsum('ji,ji->i', self.X, self.X)
        self.Bhat = (self.X.T @ self.Y) / XtX_vec[:,np.newaxis]
        # Compute se(betahat)
        Xr = self.Y - np.einsum('ij,jk->jik', self.X, self.B)
        Re = np.einsum('ijk,ijk->ik', Xr, Xr)
        self.Shat = np.sqrt(Re / XtX_vec[:,np.newaxis] / (self.X.shape[0] - 2))

    def remove_covariates(self):
        if self.Z is not None:
            self.Y -= self.Z @ (np.linalg.inv(self.Z.T @ self.Z) @ self.Z.T @ self.Y)
            self.Z = None

    def center_data(self):
        if self.X is not None and not self.x_centered:
            self.X -= np.mean(self.X, axis=0, keepdims=True)
            self.x_centered = True
        if self.Y is not None and not self.y_centered:
            self.Y -= np.mean(self.Y, axis=0, keepdims=True)
            self.y_centered = True
        if self.Z is not None and not self.z_centered:
            self.Z -= np.mean(self.Z, axis=0, keepdims=True)
            self.z_centered = True

    def set_xcorr(self, xcorr_file):
        if os.path.isfile(xcorr_file):
            self.xcorr = load_rds(xcorr_file)
        else:
            self.xcorr = np.corrcoef(self.X, rowvar = False)
            self.xcorr = (np.square(self.xcorr) * np.sign(self.xcorr)).astype(np.float16)

    def plot_xcorr(self, out):
        use_abs = np.sum((self.xcorr < 0).values.ravel()) == 0
        fig, ax = plt.subplots()
        cmap = sns.cubehelix_palette(50, hue=0.05, rot=0, light=1, dark=0, as_cmap=True)
        sns.heatmap(self.xcorr, ax = ax, cmap = cmap, vmin=-1 if not use_abs else 0,
                    vmax=1, square=True, xticklabels = False, yticklabels = False)
        ax = plt.gca()
        plt.savefig(out, dpi = 500)
        
    def permute_X_columns(self):
        '''
        Permute X columns, i.e. break blocked correlation structure
        '''
        np.random.shuffle(self.X) 
        
    def plot_property_vector(self, yaxis, zaxis, xz_cutoff = None, out = '/tmp/1.pdf',
                            conf = {'title': '', 'ylabel': '', 'zlabel': ''}):
        '''
        - yaxis can be eg $\beta$ or log10BF or -log10Prob
        - zaxis can be some other quantity whose value will be 
        reflected by color shade
        - xz_cutoff: (c1, c2). c1 is correlation cutoff to highlight
        when c2 is satisfied by a given position on x-axis
        '''
        xaxis = [x+1 for x in range(len(yaxis))]
        cmap = sns.cubehelix_palette(start=2.8, rot=.1, as_cmap=True)
        f, ax = plt.subplots(figsize=(18,5))
        points = ax.scatter(xaxis, yaxis, c=zaxis, cmap=cmap)
        f.colorbar(points, label=conf['zlabel'])
        if xz_cutoff is not None:
            c1, c2 = xz_cutoff
            for idx, item in enumerate(zaxis):
                if item > c2:
                    ax.scatter(xaxis[idx], yaxis[idx], s=80, 
                               facecolors='none', edgecolors='r')
                    for ii, xx in enumerate(self.corr[idx,:]):
                        if xx > c1 and xx < 1.0:
                            ax.scatter(xaxis[ii], yaxis[ii], 
                                       color='y', marker='+')
        ax.set_title(conf['title'])
        ax.set_ylabel(conf['ylabel'])
        plt.gca()
        plt.savefig(out, dpi = 500)
        
    def get_representative_features(self, block_r2 = 0.8, block_size = 10, max_indep_r2 = 0.02):
        '''
        Based on xcorr matrix, select "most representative features". 
        That is, these features are potentially most convoluted by other features (have stronger xcorr)
        yet are independent among each other.
        - block_r2: definition of correlated block -- abs squared correlation have to be > cutoff1
        - block_size: define a large enough block -- block size have to be > block_size
        - max_indep_r2: now select features that are completely independent -- r2 < max_indep_r2
        '''
        if self.xcorr is None:
            self.set_xcorr('')
        # get r2 summary
        r2 = pd.DataFrame(self.xcorr)
        strong_r2_count = ((np.absolute(r2) > block_r2) * r2).sum(axis = 0).sort_values(ascending = False)
        strong_r2_count = strong_r2_count[strong_r2_count > block_size]
        # filter by r2
        exclude = []
        for x in strong_r2_count.index:
            if x in exclude:
                continue
            for y in strong_r2_count.index:
                if y in exclude or y == x:
                    continue
                if np.absolute(r2[x][y]) > max_indep_r2:
                    exclude.append(y)
        return [x for x in strong_r2_count.index if not x in exclude]

    def __str__(self):
        return pformat(self.__dict__, indent = 4)
    
class UnivariateMixture:
    '''Simulated distributions of Stephens 2017 (ASH paper)'''
    def __init__(self, dim):
        self.size = dim
        self.pi0 = 0
        self.pis = []
        self.mus = []
        self.sigmas = []
        self.coef = []
        
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
        
    def get_effects(self):
        '''
        beta ~ \pi_0\delta_0 + \sum \pi_i N(0, sigma_i)
        '''
        sigmas = np.diag(self.sigmas)
        assert (len(self.pis), len(self.pis)) == sigmas.shape
        masks = np.random.multinomial(1, self.pis, size = self.size)
        mix = np.random.multivariate_normal(self.mus, sigmas, self.size)
        self.coef = np.sum(mix * masks, axis = 1) * np.random.binomial(1, 1 - self.pi0, self.size)
        
    def swap_top_effects(self, top_index):
        '''Set top effects to top indices'''
        nb = [0] * len(self.coef)
        beta = sorted(self.coef, key=abs, reverse=True)
        for idx in top_infex:
            nb[idx] = beta.pop(0)
        random.shuffle(beta)
        for idx in range(len(nb)):
            if not idx in top_index:
                nb[idx] = beta.pop(0)
        assert len(beta) == 0
        self.coef = np.array(nb)
        
    def __str__(self):
        params = ' + '.join(["{} N({}, {}^2)".format(x,y,z) for x, y, z in zip(self.pis, self.mus, self.sigmas)])
        return '{:.3f} \delta_0 + {:.3f} [{}]'.format(self.pi0, 1 - self.pi0, params)
    
class MultivariateMixture:
    '''FIXME: ideally implement Urbut 2017 simulated covs'''
    def __init__(self, dim):
        self.J, self.R = dim
        self.pis = dict([('null', 0)])
        self.mus = []
        self.Us = dict()
        self.coef = []
        self.grid = [0.1,0.5,1,2]
        self._init_canonical()

    def set_pi0(self, pi0):
        self.pis['null'] = pi0
        
    def set_grid(self, grid):
        self.grid = grid
        
    def _init_canonical(self):
        '''
        U is a dict of 
        - "identity" for the identity (effects are independent among conditions);
        - "singletons" for the set of matrices with just one non-zero entry x_{jj} = 1 (j=1,...,R); (effect specific to condition j);
        - "equal_effects" for the matrix of all 1s (effects are equal among conditions);
        - "simple_het" for a set of matrices with 1s on the diagonal and all off-diagonal elements equal to pho; (effects are correlated among conditions).
        '''
        pho = [0.25, 0.5, 0.75]
        self.Us['null'] = np.zeros((self.R, self.R))
        self.Us['identity'] = np.identity(self.R)
        for i in range(self.R):
            self.Us[f'singleton_{i+1}'] = np.diagflat([1 if idx == i else 0 for idx in range(self.R)])
        self.Us['equal_effects'] = np.ones((self.R, self.R))
        for idx, item in enumerate(sorted(pho)):
            self.Us[f'simple_het_{idx+1}'] = np.ones((self.R, self.R)) * item
            np.fill_diagonal(self.Us[f'simple_het_{idx+1}'], 1)
            
    def set_shared(self):
        '''
        All weights are on equal effects
        '''
        self.pis['equal_effects'] = 1 - self.pis['null']
        for k in self.Us:
            if k not in self.pis:
                self.pis[k] = 0
                
    def set_low_het(self):
        '''
        All weights are on small het effects
        '''
        self.pis['simple_het_1'] = 1 - self.pis['null']
        for k in self.Us:
            if k not in self.pis:
                self.pis[k] = 0
                
    def set_indep(self):
        '''
        All weights are on identity effects
        '''
        self.pis['identity'] = 1 - self.pis['null']
        for k in self.Us:
            if k not in self.pis:
                self.pis[k] = 0

    def set_singleton(self, index):
        '''
        All weights evenly set to given index of singleton effects
        '''
        index = [int(x) for x in index if x <= self.R and x > 1]
        weight = (1 - self.pis['null']) / len(index)
        for item in index:
            self.pis[f'singleton_{item}'] = weight
        for k in self.Us:
            if k not in self.pis:
                self.pis[k] = 0        
        
    def apply_grid(self):
        def product(x,y):
            for item in y:
                yield x*item
        self.Us = dict(sum([[(f"{p}_{i+1}", g) for i, g in enumerate(product(self.Us[p], np.square(self.grid)))] for p in self.Us if p != 'null'], []) + \
                      [('null', self.Us['null'])])
        nG = len(self.grid)
        for k in list(self.pis.keys()):
            if k == 'null':
                continue
            for g in range(nG):
                self.pis[f'{k}_{g+1}'] = self.pis[k] / nG
            del self.pis[k]
