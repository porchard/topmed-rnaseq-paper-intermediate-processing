#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from cyvcf2 import VCF


def get_alt_allele_counts(vcf_path, region=None, variants=None, samples=None):
    if region is not None and variants is not None:
        raise ValueError('Only region or variants can be set, not both')
    
    vcf = VCF(vcf_path, gts012=True, samples=samples)
    alt_allele_counts = []
    variants_pulled = []

    if variants is not None:
        for variant in variants:
            chrom, pos = variant.split('_')[:2]
            start, end = int(pos) - 1, int(pos)
            for v in vcf(f'{chrom}:{start}-{end}'):
                if v.ID == variant:
                    alt_allele_counts.append(list(v.gt_types))
                    variants_pulled.append(v.ID)
    elif region is not None:
        for v in vcf(f'{chrom}:{start}-{end}'):
            alt_allele_counts.append(list(v.gt_types))
            variants_pulled.append(v.ID)
    else:
        for v in vcf():
            alt_allele_counts.append(list(v.gt_types))
            variants_pulled.append(v.ID)

    alt_allele_counts_df = pd.DataFrame(alt_allele_counts, columns=vcf.samples, index=variants_pulled)
    return alt_allele_counts_df


def alt_allele_counts_matrix_to_r2_matrix(x):
    """
    Rows are variants, in format {chrom}_{pos}_{ref}_{alt}
    Columns are samples
    """

    # split by chrom; all between-chrom values are set to 0
    r2_by_chrom = []
    for chrom, df in x.groupby(lambda x: x.split('_')[0]):
        r2_by_chrom.append(df.T.corr() ** 2)
    return pd.concat(r2_by_chrom).fillna(0)