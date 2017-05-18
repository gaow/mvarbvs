#!/usr/bin/env python3
# Adopted by Gao Wang from Francois Aguet:
# https://github.com/broadinstitute/gtex-pipeline

import numpy as np
import pandas as pd
import gzip
import subprocess
import scipy.stats as stats
import argparse
import re, os

def annotate_tissue_data(data, fsample):
    '''Save data to tissue specific tables'''
    samples = pd.read_csv(fsample, dtype=str, delimiter='\t', header=0)
    sample_dict = {}
    for row in samples[['SAMPID', 'SMTSD']].values:
        if row[1] not in sample_dict:
            sample_dict[row[1]] = []
        if row[0] in data.columns:
            sample_dict[row[1]].append(row[0])
    sample = dict((re.sub("[\W\d]+", "_", k.strip()).strip('_'), v) for k, v in sample_dict.items() if len(v))
    data = {k: data.loc[:, sample[k]] for k in sample}
    return data

def write_per_tissue_data(data, output):
    if os.path.isfile(output):
        os.remove(output)
    for k in data:
        data[k].to_hdf(output, k, mode = 'a', complevel = 9, complib = 'zlib')

def get_donors_from_vcf(vcfpath):
    """
    Extract donor IDs from VCF
    """
    with gzip.open(vcfpath) as vcf:
        for line in vcf:
            if line.decode()[:2]=='##': continue
            break
    return line.decode().strip().split('\t')[9:]


def read_gct(gct_file, donor_ids):
    """
    Load GCT as DataFrame
    First col of expression data is ENCODE gene name, 2nd col is HUGO name
    ======================================================================
    A more memory friendly version:

    head = pd.read_csv(fdata, skiprows = 2, sep = '\t', nrows = 1)
    dt = {'Description': str, 'Name': str}
    dt.update({x: dtype for x in head.columns if x not in dt})
    data = pd.read_csv(fdata, compression='gzip', skiprows=2,
                       index_col=0, header=0, dtype = dt, sep='\t').drop('Description', 1)

    """
    df = pd.read_csv(gct_file, sep='\t', skiprows=2, index_col=0)
    df.drop('Description', axis=1, inplace=True)
    df.index.name = 'gene_id'
    return df[[i for i in df.columns if '-'.join(i.split('-')[:2]) in donor_ids]]


def normalize_quantiles(M, inplace=False):
    """
    Note: replicates behavior of R function normalize.quantiles from library("preprocessCore")

    Reference:
     [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003

    Adapted from https://github.com/andrewdyates/quantile_normalize
    """
    if not inplace:
        M = M.copy()

    Q = M.argsort(axis=0)
    m,n = M.shape

    # compute quantile vector
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:,i],i]
    quantiles = quantiles / n

    for i in range(n):
        # Get equivalence classes; unique values == 0
        dupes = np.zeros(m, dtype=np.int)
        for j in range(m-1):
            if M[Q[j,i],i]==M[Q[j+1,i],i]:
                dupes[j+1] = dupes[j]+1

        # Replace column with quantile ranks
        M[Q[:,i],i] = quantiles

        # Average together equivalence classes
        j = m-1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j-dupes[j]:j+1,i]
                M[idxs,i] = np.median(M[idxs,i])
                j -= 1 + dupes[j]
        assert j == -1

    if not inplace:
        return M


def inverse_quantile_normalization(M):
    """
    After quantile normalization of samples, standardize expression of each gene
    """
    R = stats.mstats.rankdata(M,axis=1)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[1]+1))
    return Q


def normalize_expression(expression_df, counts_df, expression_threshold=0.1, count_threshold=5, min_samples=10):
    """
    Genes are thresholded based on the following expression rules:
      >=min_samples with >expression_threshold expression values
      >=min_samples with >count_threshold read counts
    """
    # donor_ids = ['-'.join(i.split('-')[:2]) for i in expression_df.columns]
    donor_ids = expression_df.columns

    # expression thresholds
    mask = ((np.sum(expression_df>expression_threshold,axis=1)>=min_samples) & (np.sum(counts_df>count_threshold,axis=1)>=min_samples)).values

    # apply normalization
    M = normalize_quantiles(expression_df.loc[mask].values, inplace=False)
    R = inverse_quantile_normalization(M)

    quant_std_df = pd.DataFrame(data=R, columns=donor_ids, index=expression_df.loc[mask].index)
    quant_df = pd.DataFrame(data=M, columns=donor_ids, index=expression_df.loc[mask].index)
    return quant_std_df, quant_df


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate normalized expression BED files for eQTL analyses')
    parser.add_argument('expression_gct', help='GCT file with expression in normalized units, e.g., TPM or FPKM')
    parser.add_argument('counts_gct', help='GCT file with read counts')
    parser.add_argument('vcf', help='VCF file with donor IDs')
    parser.add_argument('attributes', help='Sample attributes, corresponding to the GTF')
    parser.add_argument('prefix', help='Prefix for output file names')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    parser.add_argument('--expression_threshold', type=np.double, default=0.1, help='Selects genes with > expression_threshold expression in at least min_samples')
    parser.add_argument('--count_threshold', type=np.int32, default=5, help='Selects genes with > count_threshold reads in at least min_samples')
    parser.add_argument('--min_samples', type=np.int32, default=10, help='Minimum number of samples that must satisfy thresholds')
    args = parser.parse_args()

    print('Generating normalized expression files ... ', end='', flush=True)
    donor_ids = get_donors_from_vcf(args.vcf)
    expression_df = read_gct(args.expression_gct, donor_ids)
    counts_df = read_gct(args.counts_gct, donor_ids)

    quant_std_df, quant_df = normalize_expression(expression_df, counts_df,
        expression_threshold=args.expression_threshold, count_threshold=args.count_threshold, min_samples=args.min_samples)
    print('Save to HDF5 format, full matrix and per tissue data ...', end='', flush=True)
    prefix = os.path.join(args.output_dir, args.prefix)
    quant_std_per_tissue = annotate_tissue_data(quant_std_df, args.attributes)
    quant_per_tissue = annotate_tissue_data(quant_df, args.attributes)

    quant_std_df.to_hdf(prefix + ".qnorm.std.flat.h5",  'GTExV7', mode = 'w', complevel = 9, complib = 'zlib')
    quant_df.to_hdf(prefix + ".qnorm.flat.h5", 'GTExV7', mode = 'w', complevel = 9, complib = 'zlib')
    write_per_tissue_data(quant_per_tissue, prefix + ".qnorm.h5")
    write_per_tissue_data(quant_std_per_tissue, prefix + ".qnorm.std.h5")
    print('done.')
