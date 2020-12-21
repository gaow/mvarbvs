
datpos = unlist(strsplit(dataset, ','))
chr=datpos[1]; start=datpos[2]; end=datpos[3];

geno.file = paste0(genotype_dir, 'bloodcells_chr',
                   chr, '.', start, '.', end)
ld.file = paste0(ld_dir, 'bloodcells_chr',
                   chr, '.', start, '.', end, '.matrix')
ldeigen.file = paste0(ldeigen_dir, 'bloodcells_chr',
                   chr, '.', start, '.', end, '.rds')
varY = readRDS(varY_file)$Y.cov


