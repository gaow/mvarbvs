import subprocess
import pandas as pd

def dap(x, y, prefix, **kwargs):
    '''
    Wen et al 2016 AJHG
    '''
    with open(f'{prefix}.data', 'w') as f:
        print(*(['pheno', 'pheno', '0'] + list(y.ravel())), file=f)
        for j in range(x.shape[1]):
            print(*(['geno', 'snp{}'.format(j), '0'] + list(x[:,j])), file=f)
    with open(f'{prefix}.grid', 'w') as f:
        print(0, 1, file=f)
    out, err = subprocess.Popen(['dap', '-d', f'{prefix}.data', '-g', f'{prefix}.grid', '-all'], stdout=subprocess.PIPE).communicate()
    with io.StringIO(str(out, 'utf-8')) as f:
        for line in f:
            if line.startswith('Posterior inclusion probability'):
                next(f)
                return pd.read_table(f, header=None, names=['rank', 'snp', 'pip', 'bf'], sep='\s+').set_index('snp')
    raise RuntimeError('Failed to parse output')
