import itertools
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
COLORS = ['#348ABD', '#7A68A6', '#A60628', '#467821', '#FF0000', '#188487', '#E2A233',
              '#A9A9A9', '#000000', '#FF00FF', '#FFD700', '#ADFF2F', '#00FFFF']
color_mapper = np.vectorize(lambda x: dict([(i,j) for i,j in enumerate(COLORS)]).get(x))

def get_set_positions(in_cs):
    positions = []
    curr_len = 0
    for b, g in itertools.groupby(in_cs):
        lg = len(list(g))
        if b:
            positions.append((curr_len, curr_len + lg))
        curr_len += lg
    return positions

def get_cs(in_cs):
    return [get_set_positions(c) for c in in_cs]

def get_purity_values(cs, ld_mat):
    lds = []
    for item in cs:
        tmp = []
        for ii in item:
            value = np.array(np.absolute(ld_mat[ii[0]:ii[1],ii[0]:ii[1]])).flatten()
            value = [x for x in value if x < 1]
            if len(value):
                tmp.append((np.min(value), np.mean(value), np.median(value), np.max(value)))
            else:
                tmp.append((1,1,1,1))
        lds.append(np.array(tmp))
    return lds, {'min': 0, 'mean': 1, 'median': 2, 'max': 3}

def plot_sets(in_cs, n_in_cs, lfsr, true_coef, ld_mat, seg_prefix):
    ld_status = dict()
    for idx, k in enumerate(in_cs.keys()):
        purity = plot_one_set(np.array(in_cs[k]), 
                        np.array(n_in_cs)[:, idx], 
                        np.array(lfsr)[:, idx], 
                        np.array(true_coef)[:, idx], 
                        ld_mat,
                        f'{seg_prefix}.{idx+1}.png')
        ld_status[idx+1] = purity
    return ld_status
        
def plot_segments(cs, nonzeros, n_in_cs, lfsr, ld_min, seg_file):
    lengths = np.concatenate([np.full(len(x), idx) for idx, x in enumerate(cs)])
    lengths = np.array([(len(lengths) - idx, x) for idx, x in enumerate(lengths)])
    sets = np.vstack(np.array(cs))
    y, c, x1, x2 = np.hstack((lengths, sets)).T
    plt.figure(figsize=(6,4))
    plt.vlines(nonzeros, 0, max(y), colors='#dcdcdc')
    for idx, value in enumerate(cs):
        keep = [j for j, item in enumerate(c) if item == idx]
        nonzero_captured = [len([s for s in value if item >= s[0] and item <= s[1]]) > 0 for item in nonzeros]
        singletons = [item for item in value if item[1] - item[0] == 1]
        label = f'{idx+1}: {int(n_in_cs[idx])} vars, {sum(nonzero_captured)}/{len(nonzero_captured)} hit, min(LD)={ld_min[idx]:.2f}, lfsr={lfsr[idx]:.2f}'
        plt.hlines(np.take(y, keep), np.take(x1, keep), np.take(x2, keep), 
                   colors = color_mapper(np.take(c, keep)), label = label)
    plt.xlabel('')
    plt.ylabel('')
    plt.yticks([])
    plt.title("95% CI sets")
#     plt.legend(loc='upper center', bbox_to_anchor=(0, -0.1),
#                 mode="expand", borderaxespad=0, ncol=1)
    plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
    plt.savefig(seg_file, dpi=500, bbox_inches='tight')                    
    plt.gca()
    
def plot_purity(cs, ld_mat, purity_file, ld_type = 'min'):
    lds, col_map = get_purity_values(cs, ld_mat)
    ld_idx = col_map[ld_type]
    idx = 0
    plt.figure(figsize=(12, 8))
    L = len(lds)
    cols = 3
    rows = L // cols + L % cols
    position = range(1,L + 1)
    for size, ld in zip(cs, lds):
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
    
def get_cs_purity(cs, ld_mat):
    lds = []
    for item in cs:
        idx = sum([list(range(x[0], x[1])) for x in item], [])
        value = np.array(np.absolute(ld_mat[idx,:][:,idx])).flatten()
        value = [x for x in value if x < 1]
        if len(value):
            lds.append((np.min(value), np.mean(value), np.median(value), np.max(value)))
        else:
            lds.append((1,1,1,1))
    return np.array(lds)

def plot_one_set(in_cs, n_in_cs, lfsr, true_coef, ld_mat, seg_file):
        cs = get_cs(in_cs)
        nonzeros = [i for i, item in enumerate(true_coef) if item != 0]
        purity = get_cs_purity(cs, ld_mat)
        plot_segments(cs, nonzeros, n_in_cs, lfsr, purity[:,0], seg_file)
        return purity
