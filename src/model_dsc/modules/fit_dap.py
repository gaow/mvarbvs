import subprocess
import pandas as pd
import numpy as np

def dap_single(x, y, prefix, args):
    '''
    Wen et al 2016 AJHG
    Here I run DAP version 2
    relevant output:
    ================
    Posterior expected model size: 0.500 (sd = 0.500)
    LogNC = -0.30685 ( Log10NC = -0.133 )
    Posterior inclusion probability

    ((1))              7492 6.68581e-05       0.000 1
    ((2))              7490 6.68581e-05       0.000 1
    ((3))              7484 6.68581e-05       0.000 1
    ((4))              7486 6.68581e-05       0.000 1
    ((5))              7481 6.68581e-05       0.000 1
    ((6))              7476 6.68581e-05       0.000 1
    ((7))              7479 6.68581e-05       0.000 1
    ((8))              7491 6.68046e-05       0.000 2
    ((9))              7483 6.68046e-05       0.000 2
    ((10))             7485 6.68046e-05       0.000 2
    ((11))             7488 6.68046e-05       0.000 2
    ((12))             7474 6.68046e-05       0.000 2
    ((13))             7475 6.68046e-05       0.000 2
    ((14))             7478 6.68046e-05       0.000 2
    ((15))             7465 6.68046e-05       0.000 2
    ((16))             7473 6.68046e-05       0.000 2
    ((17))             7470 6.68046e-05       0.000 2
    ((18))             7467 6.68046e-05       0.000 2
    ((19))             7461 6.68046e-05       0.000 2
    ((20))             7459 6.68046e-05       0.000 2
    ((21))             7482 6.67422e-05       0.000 -1
    ((22))             7489 6.67422e-05       0.000 -1
    ((23))             7487 6.67422e-05       0.000 -1
    ((24))             7477 6.67422e-05       0.000 -1
    ((25))             7480 6.67422e-05       0.000 -1
    ((26))             7463 6.67422e-05       0.000 -1
    ...
    Independent association signal clusters

         cluster         member_snp      cluster_pip      average_r2
           {1}              7            4.680e-04          0.951                 0.951   0.037
           {2}             13            8.685e-04          0.623                 0.037   0.623

    '''
    names = np.array([('geno', i+1, 0) for i in range(x.shape[1])])
    with open(f'{prefix}.data', 'w') as f:
        print(*(['pheno', 'pheno', '0'] + list(y.ravel())), file=f)
        np.savetxt(f, np.hstack((names, x.T)), fmt = '%s', delimiter = ' ')
    with open(f'{prefix}.grid', 'w') as f:
        print(0, 1, file=f)
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
            clusters.append([float(line[2]), float(line[3])])
    pips = pd.DataFrame(pips, columns = ['x', 'pip', 'bf', 'cluster'])
    clusters = pd.DataFrame(clusters, columns = ['cluster_pip', 'avg_r2'])
    return {'pips': pips, 'clusters': clusters}

def dap_batch(X, Y, prefix, *args):
    return dict([(r, dap_single(X, Y[:,r], f'{prefix}_condition_{r+1}', args)) for r in range(Y.shape[1])])
