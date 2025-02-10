#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import pybedtools as bt
from topmed_manuscript_clean import phenotype_id_to_gene_id, top_pip_variants
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')

parser = argparse.ArgumentParser()
parser.add_argument('--maf', required=True, nargs='+')
parser.add_argument('--ld', required=True, nargs='+')
parser.add_argument('--phenotypes', required=True)
parser.add_argument('--credible-sets', required=True)
parser.add_argument('--prefix', required=True)
args = parser.parse_args()


# load credible sets
# load phenotypes
# load MAF
# load LD

#import glob
#MAF = glob.glob('/net/topmed10/working/porchard/rnaseq/work/matching-statistics/freeze-2b/cis-eqtl/maf001/results/maf-and-ld/Whole_blood.chr*.maf.txt')
#LD = glob.glob('/net/topmed10/working/porchard/rnaseq/work/matching-statistics/freeze-2b/cis-eqtl/maf001/results/maf-and-ld/Whole_blood.chr*.ld.txt')
#PHENOTYPES = '/net/topmed10/working/porchard/rnaseq/work/freezes/freeze-1.1RNA/files/cis-eqtl/tensorqtl-input/Whole_blood.tensorqtl-in.phenotypes.bed.gz'
#CREDIBLE_SETS = '/net/topmed10/working/porchard/rnaseq/work/freezes/freeze-1.1RNA/files/cis-eqtl/susie/Whole_blood.maf001.cs.txt'

MAF = args.maf
LD = args.ld
PHENOTYPES = args.phenotypes
CREDIBLE_SETS = args.credible_sets
PREFIX = args.prefix


maf = pd.concat([pd.read_csv(f, delim_whitespace=True) for f in MAF])

ld = pd.concat([pd.read_csv(f, delim_whitespace=True) for f in LD])

# dict of variant --> LD buddies
ld_buddies = {variant: [] for variant in maf.SNP}
for variant_a, variant_b in zip(ld.SNP_A, ld.SNP_B):
    ld_buddies[variant_a].append(variant_b)
    ld_buddies[variant_b].append(variant_a)


# calculate the number of genes each variant was tested against
variant_df = maf[['SNP']]
variant_df[['chrom', 'pos']] = variant_df['SNP'].str.split('_', expand=True)[[0, 1]]
variant_df['start'] = variant_df.pos.astype(int) - 1
variant_df['end'] = variant_df.start + 1
variant_df = variant_df[['chrom', 'start', 'end', 'SNP']]

phenotypes = pd.read_csv(PHENOTYPES, sep='\t', usecols=[0, 1, 2, 3])
phenotypes.columns = ['chrom', 'start', 'end', 'phenotype']
phenotypes['gene'] = phenotypes.phenotype.map(phenotype_id_to_gene_id) # for splicing phenotypes, want genes, not introns
phenotypes = phenotypes[['chrom', 'start', 'end', 'gene']].drop_duplicates()


SCAN_WINDOW = 1000000
gene_bt = bt.BedTool.from_dataframe(phenotypes).sort()
variant_bt = bt.BedTool.from_dataframe(variant_df).sort()
gene_to_variants = gene_bt.window(variant_bt, w=SCAN_WINDOW).to_dataframe()
variant_to_genes_tested = gene_to_variants.thickEnd.value_counts().to_dict()


variant_stats = maf[['SNP', 'MAF']].rename(columns={'SNP': 'variant', 'MAF': 'maf'})
variant_stats['n_ld_buddies'] = variant_stats.variant.map(lambda x: len(ld_buddies[x]))
variant_stats['n_genes_tested'] = variant_stats.variant.map(lambda x: variant_to_genes_tested[x] if x in variant_to_genes_tested else 0)


# for each credible set, get the TOP pip variant
# then select a control credible set by matching on MAF, number of LD buddies, and number of genes tested against
susie = pd.read_csv(CREDIBLE_SETS, sep='\t')


susie['maf'] = np.minimum(susie.af, 1-susie.af)
susie_cs_size = susie.groupby(['phenotype_id', 'cs_id']).size().rename('n_variants').reset_index()
susie_top_pip_variant = top_pip_variants(susie)
susie_top_pip_variant = susie_top_pip_variant.merge(susie_cs_size)
susie_top_pip_variant['n_genes_tested'] = susie_top_pip_variant.variant_id.map(lambda x: variant_to_genes_tested[x] if x in variant_to_genes_tested else 0)


# this is wrong, but not yet sure why. Plots just look bad.
# control_variants = []

# # use nested loop to speed up subsetting
# # could also do by chrom
# for N_VARIANTS, susie_subset in susie_top_pip_variant.groupby('n_variants'):
#     LD_BUDDY_BUFFER = 0 if N_VARIANTS < 10 else 0.1 * N_VARIANTS
#     variant_stats_subset = variant_stats[(variant_stats.n_ld_buddies-(N_VARIANTS-1)).abs() <= LD_BUDDY_BUFFER]
#     for N_GENES_TESTED, MAF in zip(susie_subset.n_genes_tested, susie_subset.maf):
#         N_GENES_TESTED_BUFFER = max(4, 0.2 * N_GENES_TESTED)
#         MIN_N_GENES_TESTED = 1
#         MAF_BUFFER = max(MAF*0.1, 0.005) # allow more wiggle room for higher MAF variants. Alternatively, just choose the one with closest MAF match?

#         if len(control_variants) % 100 == 0:
#             logging.info('Processed {:,} control variants'.format(len(control_variants)))
#         tmp = variant_stats_subset[variant_stats_subset.n_genes_tested >= MIN_N_GENES_TESTED]
#         tmp = tmp[(tmp.n_genes_tested-N_GENES_TESTED).abs() <= N_GENES_TESTED_BUFFER]
#         tmp['maf_delta'] = (tmp.maf-MAF).abs()
#         tmp = tmp[~tmp.variant.isin(susie.variant_id)]
#         tmp = tmp[~tmp.variant.isin(control_variants)] # TODO: do we care about LD buddies?
#         if len(tmp) == 0:
#             logging.info('No control variants found for {:,} variants, {:,} genes tested, MAF {:.3f}'.format(N_VARIANTS, N_GENES_TESTED, MAF))
#             control_variants.append('')
#         else:
#             control_variants.append(tmp[tmp.maf_delta==tmp.maf_delta.min()].sample().variant.values[0])



control_variants = []


for N_VARIANTS, N_GENES_TESTED, MAF in zip(susie_top_pip_variant.n_variants, susie_top_pip_variant.n_genes_tested, susie_top_pip_variant.maf):
    
    LD_BUDDY_BUFFER = 0 if N_VARIANTS < 10 else 0.1 * N_VARIANTS
    N_GENES_TESTED_BUFFER = max(4, 0.2 * N_GENES_TESTED)
    MIN_N_GENES_TESTED = 1
    MAF_BUFFER = max(MAF*0.1, 0.005) # allow more wiggle room for higher MAF variants. Alternatively, just choose the one with closest MAF match?

    if len(control_variants) % 100 == 0:
        logging.info('Processed {:,} control variants'.format(len(control_variants)))
    tmp = variant_stats[(variant_stats.n_ld_buddies-(N_VARIANTS-1)).abs() <= LD_BUDDY_BUFFER]
    tmp = tmp[tmp.n_genes_tested >= MIN_N_GENES_TESTED]
    tmp = tmp[(tmp.n_genes_tested-N_GENES_TESTED).abs() <= N_GENES_TESTED_BUFFER]
    tmp['maf_delta'] = (tmp.maf-MAF).abs()
    tmp = tmp[~tmp.variant.isin(susie.variant_id)]
    tmp = tmp[~tmp.variant.isin(control_variants)] # TODO: do we care about LD buddies?
    if len(tmp) == 0:
        logging.info('No control variants found for {:,} variants, {:,} genes tested, MAF {:.3f}'.format(N_VARIANTS, N_GENES_TESTED, MAF))
        control_variants.append('')
    else:
        control_variants.append(tmp[tmp.maf_delta==tmp.maf_delta.min()].sample().variant.values[0])



susie_top_pip_variant['control_variant'] = control_variants 
susie_top_pip_variant = susie_top_pip_variant.merge(variant_stats.rename(columns=lambda x: 'control_' + x), how='left')
susie_top_pip_variant.to_csv(f'{PREFIX}control-top-pip-variants.tsv', sep='\t', index=False)

control_cs = []
tmp = susie_top_pip_variant[susie_top_pip_variant.control_variant!='']
for phenotype_id, cs_id, v in zip(tmp.phenotype_id, tmp.cs_id, tmp.control_variant):
    variant_id = [v] + ld_buddies[v]
    maf = variant_stats.set_index('variant').loc[variant_id].maf.to_list()
    pip = [1/len(variant_id)] * len(variant_id)
    df = pd.DataFrame({'variant_id': variant_id, 'maf': maf, 'pip': pip, 'cs_id': cs_id, 'phenotype_id': phenotype_id})
    df = df[['phenotype_id', 'variant_id', 'pip', 'maf', 'cs_id']]
    control_cs.append(df)
control_cs = pd.concat(control_cs)
control_cs.to_csv(f'{PREFIX}control-credible-sets.tsv', sep='\t', index=False)


susie_top_pip_variant['log10_maf'] = np.log10(susie_top_pip_variant.maf)
susie_top_pip_variant['log10_control_maf'] = np.log10(susie_top_pip_variant.control_maf)
susie_top_pip_variant['control_n_variants'] = susie_top_pip_variant.control_n_ld_buddies+1


fig, axs = plt.subplots(figsize=(3*4, 4), ncols=3)

ax = axs[0]
sns.scatterplot(x='n_genes_tested', y='control_n_genes_tested', data=susie_top_pip_variant, ax=ax, alpha=0.1)

ax = axs[1]
sns.scatterplot(x='log10_maf', y='log10_control_maf', data=susie_top_pip_variant, ax=ax, alpha=0.1)

ax = axs[2]
ax.set_xscale('log')
ax.set_yscale('log')
sns.scatterplot(x='n_variants', y='control_n_variants', data=susie_top_pip_variant, ax=ax, alpha=0.1)

fig.tight_layout()
fig.savefig(f'{PREFIX}control-credible-sets.png', bbox_inches='tight', dpi=300, facecolor='white')