#!/usr/bin/env python
# coding: utf-8

import json
import sys
import os
import re
import argparse
import logging

import pandas as pd

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

def phenotype_id_to_gene_id(x):
    # first, try to match with a version
    ENSEMBL_RE_WITH_VERSION = 'ENSG\d+\.\d+'
    ENSEMBL_RE_WITHOUT_VERSION = 'ENSG\d+'
    with_version = re.search(ENSEMBL_RE_WITH_VERSION, x)
    without_version = re.search(ENSEMBL_RE_WITHOUT_VERSION, x)
    if with_version:
        return with_version.group(0)
    elif without_version:
        return without_version.group(0)
    else:
        raise ValueError(f'Not able to infer gene ID from {x}')


parser = argparse.ArgumentParser()
parser.add_argument('--colocs-only', default=False, action='store_true')
parser.add_argument('--coloc-out', required=True)
parser.add_argument('--panukbb-root', required=True)
parser.add_argument('--susie-json', nargs='+', required=True)
args = parser.parse_args()

COLOC_OUT = args.coloc_out
SUSIE_JSON = args.susie_json
PANUKBB_ROOT = args.panukbb_root


coloc = pd.read_csv(COLOC_OUT, sep='\t')

susie_info = dict()
for susie_json in SUSIE_JSON:
    with open(susie_json, 'r') as fh:
        tmp = json.load(fh)
        for k in tmp:
            if k in susie_info:
                raise ValueError(f'Duplicate key {k} in {susie_json}')
        susie_info.update(tmp)

for xqtl_info in ['modality', 'tissue', 'maf', 'phenotype', 'ancestry']:
    coloc[f'xqtl_{xqtl_info}'] = coloc.xqtl_file.map(lambda x: susie_info[x][xqtl_info])

coloc['xqtl_gene'] = coloc.xqtl_phenotype.map(phenotype_id_to_gene_id)

for gwas_info in ['region', 'ancestry', 'trait']:
    coloc[f'gwas_{gwas_info}'] = coloc.gwas_file.map(lambda x: susie_info[x][gwas_info])

coloc['gwas_signal'] = coloc.gwas_trait + '___' + coloc.gwas_ancestry + '___' + coloc.gwas_region

# A few SuSiE GWAS runs did not converge. Exclude these.
gwas_susie_convergence = pd.read_csv(os.path.join(PANUKBB_ROOT, 'work/ancestry-specific-finemapping/lift-susie/results/susie-cs-and-convergence/gwas.converged.txt'), sep='\t', header=None, names=['trait_ancestry_region', 'converged'])
gwas_susie_no_converge = gwas_susie_convergence[~gwas_susie_convergence.converged].trait_ancestry_region.unique()

coloc = coloc[~coloc.gwas_signal.isin(gwas_susie_no_converge)]

assert(all(coloc.gwas_cs=='L'+coloc.idx2.astype(str)))


def handle_multicoloc(coloc, PP4_threshold=0.8):
    # handle cases of multi-coloc, meaning:
    # cases where a single xQTL CS coloced with more than one GWAS CS for the same trait+ancestry+region
    # cases where a single GWAS CS coloced with more than one xQTL CS for the same modality+tissue+maf+ancestry+gene (not phenotype)
    # in such cases, keep the coloc with the highest posterior probability
    
    # first eliminate cases where a single GWAS CS coloced with more than one xQTL CS for the same modality+tissue+maf+ancestry+gene (not phenotype)
    logging.info('{:,} pairs tested for coloc'.format(len(coloc)))
    coloced = coloc[coloc['PP.H4.abf']>=PP4_threshold]

    logging.info('{:,} pairs coloc @ PP4 threshold = {}'.format(len(coloced), PP4_threshold))
    logging.info('Filtering cases where a single xQTL CS coloced with more than one GWAS CS for the same trait+ancestry+region or a single GWAS CS coloced with more than one xQTL CS for the same modality+tissue+maf+ancestry+gene')

    #print(coloced.groupby(['gwas_signal', 'gwas_cs', 'xqtl_modality', 'xqtl_tissue', 'xqtl_maf', 'xqtl_ancestry', 'xqtl_gene']).size().value_counts())
    coloced = coloced.sort_values(['PP.H4.abf'], ascending=False).groupby(['gwas_signal', 'gwas_cs', 'xqtl_modality', 'xqtl_tissue', 'xqtl_maf', 'xqtl_ancestry', 'xqtl_gene']).head(1)
    #print(coloced.groupby(['gwas_signal', 'gwas_cs', 'xqtl_modality', 'xqtl_tissue', 'xqtl_maf', 'xqtl_ancestry', 'xqtl_gene']).size().value_counts())

    # second, eliminate cases where a single xQTL CS coloced with more than one GWAS CS for the same trait+ancestry+region
    #print(coloced.groupby(['gwas_signal', 'xqtl_modality', 'xqtl_tissue', 'xqtl_maf', 'xqtl_ancestry', 'xqtl_phenotype', 'xqtl_cs']).size().value_counts())
    coloced = coloced.sort_values(['PP.H4.abf'], ascending=False).groupby(['gwas_signal', 'xqtl_modality', 'xqtl_tissue', 'xqtl_maf', 'xqtl_ancestry', 'xqtl_phenotype', 'xqtl_cs']).head(1)
    #print(coloced.groupby(['gwas_signal', 'xqtl_modality', 'xqtl_tissue', 'xqtl_maf', 'xqtl_ancestry', 'xqtl_phenotype', 'xqtl_cs']).size().value_counts())

    logging.info('After filtering, {:,} coloced pairs remain'.format(len(coloced)))

    return coloced


if args.colocs_only:
    coloced = handle_multicoloc(coloc)
    coloced.to_csv(sys.stdout, sep='\t', index=False)
else:
    coloc.to_csv(sys.stdout, sep='\t', index=False)