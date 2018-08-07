import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COLORS = ['#348ABD', '#7A68A6', '#A60628', '#467821', '#FF0000', '#188487', '#E2A233',
              '#A9A9A9', '#000000', '#FF00FF', '#FFD700', '#ADFF2F', '#00FFFF']
color_mapper = np.vectorize(lambda x: dict([(i,j) for i,j in enumerate(COLORS)]).get(x))

class SusieReporter:
    def __init__(self, in_cs, lfsr, true_coef, ld_mat):
        self.nonzeros = [i for i, item in enumerate(true_coef) if item != 0]
        self.ld = ld_mat        
        self.cs = self.get_cs(in_cs)
        self.lfsr = lfsr
        self.signal_captured = []
        if len(self.nonzeros):
            for value in self.cs:
                self.signal_captured.append([len([s for s in value if item >= s[0] and item <= s[1]]) > 0 for item in self.nonzeros])
        self.n_in_cs = np.sum(in_cs, axis=1)
        self.purity, self.purity_label = self.get_cs_purity()

        
    @staticmethod
    def get_set_positions(in_cs):
        positions = []
        curr_len = 0
        for b, g in itertools.groupby(in_cs):
            lg = len(list(g))
            if b:
                positions.append((curr_len, curr_len + lg))
            curr_len += lg
        return positions

    def get_cs(self, in_cs):
        return [self.get_set_positions(c) for c in in_cs]

    def get_purity_values(self):
        purity = []
        for item in self.cs:
            tmp = []
            for ii in item:
                value = np.array(np.absolute(self.ld[ii[0]:ii[1],ii[0]:ii[1]])).flatten()
                value = [x for x in value if x < 1]
                if len(value):
                    tmp.append((np.min(value), np.mean(value), np.median(value), np.max(value)))
                else:
                    tmp.append((1,1,1,1))
            purity.append(np.array(tmp))
        return purity, {'min': 0, 'mean': 1, 'median': 2, 'max': 3}
        
    def plot_segments(self, seg_file, ld_type = 'min'):
        lengths = np.concatenate([np.full(len(x), idx) for idx, x in enumerate(self.cs)])
        lengths = np.array([(len(lengths) - idx, x) for idx, x in enumerate(lengths)])
        sets = np.vstack(np.array(self.cs))
        y, c, x1, x2 = np.hstack((lengths, sets)).T
        plt.figure(figsize=(6,4))
        plt.vlines(self.nonzeros, 0, max(y), colors='#dcdcdc')
        for idx, value in enumerate(self.cs):
            keep = [j for j, item in enumerate(c) if item == idx]
            singletons = [item for item in value if item[1] - item[0] == 1]
            label = f'{idx+1}: {int(self.n_in_cs[idx])} vars, {sum(self.signal_captured[idx])}/{len(self.signal_captured[idx])} hit, min(LD)={self.purity[:,self.purity_label[ld_type]][idx]:.2f}, lfsr={self.lfsr[idx]:.2f}'
            plt.hlines(np.take(y, keep), np.take(x1, keep), np.take(x2, keep), 
                       colors = color_mapper(np.take(c, keep)), label = label)
        plt.xlabel('')
        plt.ylabel('')
        plt.yticks([])
        plt.title("95% CS")
    #     plt.legend(loc='upper center', bbox_to_anchor=(0, -0.1),
    #                 mode="expand", borderaxespad=0, ncol=1)
        plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
        plt.savefig(seg_file, dpi=500, bbox_inches='tight')                    
        plt.gca()
    
    def plot_purity(self, purity_file, ld_type = 'min'):
        lds, col_map = get_purity_values(self.cs, self.ld)
        ld_idx = col_map[ld_type]
        idx = 0
        plt.figure(figsize=(12, 8))
        L = len(lds)
        cols = 3
        rows = L // cols + L % cols
        position = range(1,L + 1)
        for size, ld in zip(self.cs, lds):
            x = np.array([s[1] - s[0] for s in size])
            y = ld[:,ld_idx]
            z = np.full(len(x), idx)
            plt.subplot(rows,cols,position[idx])
            idx += 1
            label = f'L={idx}'
            plt.scatter(x,y,c = color_mapper(z), label = label, marker = 'x')
            plt.legend()
        plt.subplots_adjust(hspace=0.3, wspace = 0.3)
        plt.suptitle(f"95% CI set sizes vs {ld_type}(abs(LD))")
        plt.savefig(purity_file, dpi=500, bbox_inches='tight')                    
        plt.gca()
    
    def get_cs_purity(self):
        lds = []
        for item in self.cs:
            idx = sum([list(range(x[0], x[1])) for x in item], [])
            value = np.array(np.absolute(self.ld[idx,:][:,idx])).flatten()
            value = [x for x in value if x < 1]
            if len(value):
                lds.append((np.min(value), np.mean(value), np.median(value), np.max(value)))
            else:
                lds.append((1,1,1,1))
        return np.array(lds), {'min': 0, 'mean': 1, 'median': 2, 'max': 3}
    
def plot_sets(in_cs, lfsr, true_coef, ld_mat, seg_prefix, save_plot):
    ld_status = dict()
    signal_status = dict()
    for idx, k in enumerate(in_cs.keys()):
        reporter = SusieReporter(np.array(in_cs[k]), 
                                 np.array(lfsr)[:, idx], 
                                 np.array(true_coef)[:, idx], 
                                 ld_mat)
        if save_plot:
            reporter.plot_segments(f'{seg_prefix}.{idx+1}.png')
        ld_status[idx+1] = reporter.purity
        signal_status[idx+1] = reporter.signal_captured
    return ld_status, signal_status
