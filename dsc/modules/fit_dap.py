import subprocess
import pandas as pd
import numpy as np

def dap_single(x, y, prefix, r, args):
    names = np.array([('geno', i+1, str(r)) for i in range(x.shape[1])])
    with open(f'{prefix}.data', 'w') as f:
        print(*(['pheno', 'pheno', str(r)] + list(y.ravel())), file=f)
        np.savetxt(f, np.hstack((names, x.T)), fmt = '%s', delimiter = ' ')
    grid = '''         
        0.0000  0.1000
        0.0000  0.2000
        0.0000  0.4000
        0.0000  0.8000
        0.0000  1.6000
        '''
    grid = '\n'.join([x.strip() for x in grid.strip().split('\n')])
    with open(f'{prefix}.grid', 'w') as f:
        print(grid, file=f)
    cmd = ['dap-g', '-d', f'{prefix}.data', '-g', f'{prefix}.grid', '-o', f'{prefix}.result', '--all'] + ' '.join(args).split()
    subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    out = [x.strip().split() for x in open(f'{prefix}.result').readlines()]
    pips = []
    clusters = []
    still_pip = True
    for line in out:
        if len(line) == 0:
            continue
        if len(line) > 2 and line[2] == 'cluster_pip':
            still_pip = False
            continue
        if still_pip and (not line[0].startswith('((') or int(line[-1]) < 0):
            continue
        if still_pip:
            pips.append([line[1], float(line[2]), float(line[3]), int(line[4])])
        else:
            clusters.append([len(clusters) + 1, float(line[2]), float(line[3])])
    pips = pd.DataFrame(pips, columns = ['snp', 'snp_prob', 'snp_log10bf', 'cluster'])
    clusters = pd.DataFrame(clusters, columns = ['cluster', 'cluster_prob', 'cluster_avg_r2'])
    clusters = pd.merge(clusters, pips.groupby(['cluster'])['snp'].apply(','.join).reset_index(), on = 'cluster')
    return {'snp': pips, 'set': clusters}

def dap_batch(X, Y, prefix, *args):
    return dict([(r, dap_single(X, Y[:,r], f'{prefix}_condition_{r+1}', r+1, args)) for r in range(Y.shape[1])])
