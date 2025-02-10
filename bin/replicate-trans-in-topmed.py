#!/usr/bin/env python
# coding: utf-8

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
import qtl.norm as norm
import statsmodels.api as sm
import topmed_genotypes as tg
import logging
import sys
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')

def residualize(x, covariates):
    assert(isinstance(x, pd.Series))
    assert(isinstance(covariates, pd.DataFrame))
    model = sm.OLS(x, sm.add_constant(covariates.T)).fit()
    return model.resid


def associate(phenotypes, genotypes, covariates):
    assert(isinstance(phenotypes, pd.Series))
    assert(isinstance(genotypes, pd.Series))
    assert(isinstance(covariates, pd.DataFrame))
    # residualize
    ADD_CONSTANT = True
    p_residuals = residualize(phenotypes, covariates)
    g_residuals = residualize(genotypes, covariates)
    p_residuals = norm.inverse_normal_transform(p_residuals)
    model = sm.OLS(p_residuals, sm.add_constant(g_residuals)) if ADD_CONSTANT else sm.OLS(p_residuals, g_residuals)
    number_samples = len(phenotypes)
    number_covariates = len(covariates)
    dof = number_samples - number_covariates - 2
    model.df_model = dof
    model.df_resid = dof
    results = model.fit()
    return results


#PHENOTYPES = '/net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/tensorqtl-in/cis-eqtl/Whole_blood.tensorqtl-in.phenotypes.bed.gz'
#COVARIATES = '/net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/tensorqtl-in/cis-eqtl/Whole_blood.tensorqtl-in.100.covariates.tsv'
#VCF_DIR = '/net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/genotypes/vcfs-updated-ids-pass-filter'
PHENOTYPES, COVARIATES, VCF_DIR, PAIRS = sys.argv[1:]

# TISSUE = os.path.basename(PHENOTYPES).split('.')[0]
# assert(TISSUE in ['Whole_blood', 'Lung', 'Nasal_epithelial', 'Monocyte', 'T_cell', 'PBMC'])

logging.info('Loading phenotypes')
phenotypes = pd.read_csv(PHENOTYPES, sep='\t').drop(columns=['#chr', 'start', 'end']).set_index('gene_id')
phenotypes.index = phenotypes.index.to_series().str.split('.', expand=True)[0] # strip gene version

logging.info('Loading covariates')
covariates = pd.read_csv(COVARIATES, sep='\t', index_col=0)
# if TISSUE == 'Whole_blood':
#     covariates = covariates.drop(index=[f'phenotype_PC{i}' for i in [25] + list(range(51, 101))])

assert(all(covariates.columns == phenotypes.columns))


pairs = pd.read_csv(PAIRS, sep='\t', header=None, names=['variant_id', 'gene_id'])
if '.' in pairs.gene_id.values[0]:
    pairs.gene_id = pairs.gene_id.str.split('.', expand=True)[0] # strip gene version
all_variants = set(pairs.variant_id)
all_variants_df = pd.DataFrame(list(all_variants), columns=['variant_id'])
all_variants_df[['chrom', 'pos']] = all_variants_df.variant_id.str.split('_', expand=True)[[0, 1]]


logging.info('Fetching {:,} variants'.format(len(all_variants_df)))
aac = pd.concat(
    [tg.get_alt_allele_counts(vcf_path=f'{VCF_DIR}/{chrom}.vcf.gz', variants=df.variant_id.unique()) for chrom, df in all_variants_df.groupby('chrom')]
).astype(int)
aac = aac[covariates.columns.to_list()]


results = []
skip_count = 0
line_count = 0

for variant, gene in zip(pairs.variant_id, pairs.gene_id):
    line_count += 1
    if line_count % 1000 == 0:
        logging.info('Processed {:,} of {:,} pairs'.format(line_count, len(pairs)))
    if gene not in phenotypes.index or variant not in aac.index:
        logging.info(f'Skipping pair {gene} - {variant}; the gene or variant is missing from input data')
        skip_count += 1
        continue
    mod = associate(phenotypes.loc[gene], aac.loc[variant], covariates)
    slope = mod.params[0]
    p = mod.pvalues[0]
    results.append([variant, gene, slope, p])

results = pd.DataFrame(results, columns=['variant_id', 'gene_id', 'slope', 'p'])
results.to_csv(sys.stdout, sep='\t', index=False)

logging.info('Done. Skipped {:,} pairs.'.format(skip_count))