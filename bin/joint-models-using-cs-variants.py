#!/usr/bin/env python
# coding: utf-8


import os
import pandas as pd
from topmed_manuscript_clean import get_alt_allele_counts, top_pip_variants
import qtl.norm as norm
import statsmodels.api as sm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--phenotypes', required=True)
parser.add_argument('--covariates', required=True)
parser.add_argument('--vcf-path', required=True)
parser.add_argument('--credible-sets', required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()

#PHENOTYPES = '/net/topmed10/working/porchard/rnaseq/work/joint-model-with-cs-variants/freeze-2b/data/phenotypes/Monocyte.ciseqtl.phenotypes.bed.gz'
#COVARIATES = '/net/topmed10/working/porchard/rnaseq/work/joint-model-with-cs-variants/freeze-2b/data/covariates/Monocyte.ciseqtl.covariates.tsv'
##VCF_PATH = '/net/topmed10/working/porchard/rnaseq/work/subset-topmed-bcf/freeze-alpha/results/vcfs-updated-ids-pass-filter'
#CREDIBLE_SETS = '/net/topmed10/working/porchard/rnaseq/work/joint-model-with-cs-variants/freeze-2b/data/cs/Monocyte.maf001.ciseqtl.cs.txt'
#SAMPLES = '/net/topmed10/working/porchard/rnaseq/work/joint-model-with-cs-variants/freeze-2b/data/samples/Monocyte.samples.txt'
#VCF_PATH = '/net/topmed10/working/porchard/rnaseq/work/cs-variant-vcf-files/freeze-2b/results/vcfs-by-chrom'
#OUTPUT = 'TEST.tsv'

PHENOTYPES = args.phenotypes
COVARIATES = args.covariates
VCF_PATH = args.vcf_path
CREDIBLE_SETS = args.credible_sets
OUTPUT = args.output


def residualize(x, covariates):
    assert(isinstance(x, pd.Series) or isinstance(x, pd.DataFrame))
    assert(isinstance(covariates, pd.DataFrame))
    if isinstance(x, pd.Series):
        model = sm.OLS(x, sm.add_constant(covariates.T)).fit()
        return model.resid
    else:
        return x.transform(lambda y: residualize(y, covariates))
    


def associate(phenotypes, genotypes, covariates):
    assert(isinstance(phenotypes, pd.Series))
    assert(isinstance(genotypes, pd.Series) or isinstance(genotypes, pd.DataFrame))
    assert(isinstance(covariates, pd.DataFrame))
    # residualize
    p_residuals = residualize(phenotypes, covariates)
    g_residuals = residualize(genotypes, covariates)
    p_residuals = norm.inverse_normal_transform(p_residuals)
    model = sm.OLS(p_residuals, sm.add_constant(g_residuals))
    number_samples = len(phenotypes)
    number_covariates = len(covariates)
    dof = number_samples - number_covariates - 2 if isinstance(genotypes, pd.Series) else number_samples - number_covariates - 1 - len(genotypes.columns) # TODO: double check this
    model.df_model = dof
    model.df_resid = dof
    results = model.fit()
    return results



cs = pd.read_csv(CREDIBLE_SETS, sep='\t')
cs = top_pip_variants(cs)


# get the genotypes
variants_to_fetch = cs[['variant_id']].drop_duplicates()
variants_to_fetch['chrom'] = variants_to_fetch['variant_id'].str.split('_', expand=True)[0]
alt_allele_counts = pd.concat([get_alt_allele_counts(os.path.join(VCF_PATH, f'{chrom}.vcf.gz'), variants=df.variant_id.unique().tolist()) for chrom, df in variants_to_fetch.groupby('chrom')])

phenotypes = pd.read_csv(PHENOTYPES, sep='\t', index_col=3).drop(columns=['#chr', 'start', 'end'])
covariates = pd.read_csv(COVARIATES, sep='\t', index_col=0)
samples = covariates.columns.to_list()
alt_allele_counts = alt_allele_counts[samples]
alt_allele_counts = alt_allele_counts.T.astype(int)








results = []
for phenotype, phenotype_data in cs.groupby('phenotype_id'):
    variants = phenotype_data.variant_id.unique().tolist()
    model_results = associate(phenotypes.loc[phenotype,samples], alt_allele_counts[variants], covariates[samples])
    results_df = pd.DataFrame(model_results.summary2().tables[1].iloc[1:, :]).reset_index().rename(columns={'index': 'variant_id'}).assign(phenotype_id=phenotype)
    results.append(results_df)
results = pd.concat(results)


# save results to tsv file
results.to_csv(OUTPUT, sep='\t', index=False)

# TODO: DF model is different (equal to # covariates + # variants??) but DF residuals is correct
# sm.OLS(phenotypes.loc[GENE,samples], sm.add_constant(pd.concat([alt_allele_counts[variants], covariates[samples].T], axis=1))).fit().summary()

