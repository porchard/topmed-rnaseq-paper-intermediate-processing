#!/usr/bin/env python3.8
# coding: utf-8


import pandas as pd
from tensorqtl import post
from qtl import norm
from pandas_plink import read_plink
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--tor', help='List of TORs for the samples in the eQTL scan')
parser.add_argument('--rna-counts', dest='counts', help='RNA count matrix')
parser.add_argument('--metadata')
parser.add_argument('--covariates')
parser.add_argument('--credible-sets', dest='cs', default=None, help='SuSiE credible set file')
parser.add_argument('--gene-variant-pairs', dest='gene_variant_pairs', default=None, help='File of gene ID, variant ID (header present); can be given instead of --credible-sets')
parser.add_argument('--plink-beds', dest='plink', nargs='+')
parser.add_argument('--out')
args = parser.parse_args()

TOR_LIST = args.tor
RNA_COUNTS = args.counts
METADATA = args.metadata
COVARIATES = args.covariates
CREDIBLE_SETS = args.cs
GENE_VARIANT_PAIRS = args.gene_variant_pairs
PLINK_BEDS = args.plink
OUT = args.out

assert(GENE_VARIANT_PAIRS is None or CREDIBLE_SETS is None)
assert(GENE_VARIANT_PAIRS is not None or CREDIBLE_SETS is not None)

CHROM_TO_PLINK_BED = {os.path.basename(f).split('.')[0]: f for f in PLINK_BEDS}

def get_genotypes(plink_bed, snps):
    (bim, fam, bed) = read_plink(plink_bed, verbose=False)
    bed = 2 - bed  # get ALT allele dosages; PLINK uses REF
    bim = bim[bim.snp.isin(set(snps))]
    snp_indices = bim.index.to_list()
    df = pd.DataFrame(bed[snp_indices,:].compute(), index=bim.snp, columns=fam.iid)
    return df


def top_pip_variants(cs):
    """
    Given a dataframe representing CS (having columns ['phenotype_id', 'variant_id', 'cs_id', 'pip']; can have other columns too),
    return a dataframe containing only the top PIP variant per credible set (credible set defined by phenotype_id + csID pair)
    """
    # Validate input
    if not isinstance(cs, pd.DataFrame):
        raise TypeError('cs must be a DataFrame')
    for i in ['phenotype_id', 'variant_id', 'cs_id', 'pip']:
        if not i in cs.columns.to_list():
            raise ValueError(f'cs must include column {i}')

    return cs.sort_values(['pip', 'variant_id'], ascending=[False, True]).groupby(['phenotype_id', 'cs_id']).head(1)


metadata = pd.read_csv(METADATA, sep='\t')
TOR_TO_NWD = dict(zip(metadata.tor, metadata.wgs))


# determine samples to use
samples = pd.read_csv(TOR_LIST, sep='\t', header=None)[0].to_list()

# get normalized counts for these samples
counts = pd.read_hdf(RNA_COUNTS)
counts = counts.loc[samples,:]


# deseq normalization
normalized_counts = norm.deseq2_normalized_counts(counts.T)
normalized_counts.columns = [TOR_TO_NWD[i] for i in normalized_counts.columns]


covariates = pd.read_csv(COVARIATES, sep='\t', index_col=0)

if CREDIBLE_SETS:
    cs = pd.read_csv(CREDIBLE_SETS, sep='\t')
    top_pip_snps = top_pip_variants(cs)
    top_pip_snps['gene_id'] = top_pip_snps.phenotype_id
else:
    top_pip_snps = pd.read_csv(GENE_VARIANT_PAIRS, sep='\t')
    top_pip_snps.columns = ['gene_id', 'variant_id']


fetch_variants = top_pip_snps.variant_id.str.split('_', expand=True)
fetch_variants.columns = ['chrom', 'pos', 'ref', 'alt']
fetch_variants['variant_id'] = top_pip_snps.variant_id


genotypes = pd.concat([get_genotypes(CHROM_TO_PLINK_BED[chrom], df.variant_id.unique()) for chrom, df in fetch_variants.groupby('chrom')])
genotypes = genotypes[normalized_counts.columns.to_list()]
variants_df = pd.DataFrame({'variant_id': genotypes.index.to_list()})
variants_df[['chrom', 'pos']] = variants_df.variant_id.str.split('_', expand=True)[[0, 1]]
variants_df.pos = variants_df.pos.astype(int)
variants_df = variants_df.set_index('variant_id')[['chrom', 'pos']]


assert(genotypes.columns == normalized_counts.columns)
assert(genotypes.columns == covariates.columns)


afc = post.calculate_afc(top_pip_snps[['gene_id', 'variant_id']].drop_duplicates(), normalized_counts, genotypes, variants_df, covariates.T, select_covariates=False)
afc.to_csv(OUT, sep='\t', index=False)