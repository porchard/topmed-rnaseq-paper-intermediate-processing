#!/usr/bin/env python
# coding: utf-8


import os
import pandas as pd
import numpy as np
import glob
import topmed_manuscript_clean as tm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('top_eqtl_glob', type=str)
parser.add_argument('top_sqtl_glob', type=str)
parser.add_argument('vcf_path', type=str)
parser.add_argument('metadata', type=str)
args = parser.parse_args()

TOP_EQTL_GLOB = args.top_eqtl_glob
TOP_SQTL_GLOB = args.top_sqtl_glob
VCF_PATH = args.vcf_path
metadata = pd.read_csv(args.metadata, sep='\t')
PREFIX = 'clump-trans-signals.'


def genotype_matrix_to_alt_allele_counts(x):
    d = {'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2}
    assert(all(np.isin(x.values.flatten(), list(d.keys()))))
    return x.replace(d)
    


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


def clump(variants, ld, ld_threshold=0.8):
    # variant must be sorted, e.g. so that the most significant one comes first
    variant_to_new_variant = {}
    for variant in variants:
        if variant not in variant_to_new_variant:
            # get variants that are in high LD and that haven't been re-assigned already
            high_ld_variants = {other_variant: pairwise_ld for other_variant, pairwise_ld in ld[variant].items() if (pairwise_ld>=ld_threshold and other_variant not in variant_to_new_variant)}
            for k in high_ld_variants:
                variant_to_new_variant[k] = variant
        else:
            continue
    variant_to_new_variant = {k: v for k, v in variant_to_new_variant.items() if k in set(variants)}
    return variant_to_new_variant



TOP_EQTL = glob.glob(TOP_EQTL_GLOB)
TOP_SQTL = glob.glob(TOP_SQTL_GLOB)

trans_eqtl = pd.concat([pd.read_csv(f, sep='\t').assign(tissue=os.path.basename(f).split('.')[0]) for f in TOP_EQTL])
trans_eqtl = trans_eqtl[trans_eqtl.qvalue<=0.05]

trans_sqtl = pd.concat([pd.read_csv(f, sep='\t').assign(tissue=os.path.basename(f).split('.')[0]) for f in TOP_SQTL])
trans_sqtl = trans_sqtl[trans_sqtl.qvalue<=0.05]

# get genotypes (alt allele counts) for all variants
all_variants = list(set(trans_eqtl.variant_id.to_list() + trans_sqtl.variant_id.to_list()))
all_variants = pd.DataFrame({'variant_id': all_variants})
all_variants['chrom'] = all_variants.variant_id.str.split('_', expand=True)[0]
genotypes = pd.concat([tm.get_genotypes_from_vcf(f'{VCF_PATH}/{chrom}.vcf.gz', variants=df.variant_id.to_list()) for chrom, df in all_variants.groupby('chrom')])
alt_allele_counts = genotype_matrix_to_alt_allele_counts(genotypes)


clump_dict = {} # modality --> tissue --> clump

for modality in ['trans_eqtl', 'trans_sqtl']:
    trans_df = trans_eqtl if modality == 'trans_eqtl' else trans_sqtl
    clump_dict[modality] = dict()
    for tissue, tissue_trans_df in trans_df.groupby('tissue'):
        tissue_metadata = metadata[(metadata.used_for_scan) & (metadata.tissue==tissue)]
        tissue_genotypes = alt_allele_counts[tissue_metadata.wgs.to_list()]
        tissue_r2 = alt_allele_counts_matrix_to_r2_matrix(tissue_genotypes).to_dict()
        clump_dict[modality][tissue] = clump(tissue_trans_df.sort_values('pval_beta', ascending=True).variant_id.unique(), tissue_r2, ld_threshold=0.9)

trans_eqtl['clumped_variant_id'] = [clump_dict['trans_eqtl'][tissue][variant] for tissue, variant in zip(trans_eqtl.tissue, trans_eqtl.variant_id)]
trans_eqtl.to_csv(f'{PREFIX}significant-trans-eqtl-clumped.tsv', sep='\t', index=False)

trans_sqtl['clumped_variant_id'] = [clump_dict['trans_sqtl'][tissue][variant] for tissue, variant in zip(trans_sqtl.tissue, trans_sqtl.variant_id)]
trans_sqtl.to_csv(f'{PREFIX}significant-trans-sqtl-clumped.tsv', sep='\t', index=False)