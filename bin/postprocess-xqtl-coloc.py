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
parser.add_argument('--susie-json', nargs='+', required=True)
args = parser.parse_args()

COLOC_OUT = args.coloc_out
SUSIE_JSON = args.susie_json



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
    coloc[f'xqtl1_{xqtl_info}'] = coloc.xqtl1_file.map(lambda x: susie_info[x][xqtl_info])
    coloc[f'xqtl2_{xqtl_info}'] = coloc.xqtl2_file.map(lambda x: susie_info[x][xqtl_info])
coloc['xqtl1_gene'] = coloc.xqtl1_phenotype.map(phenotype_id_to_gene_id)
coloc['xqtl2_gene'] = coloc.xqtl2_phenotype.map(phenotype_id_to_gene_id)

assert(all(coloc.xqtl1_tissue == coloc.xqtl2_tissue))



def handle_multicoloc(coloc, PP4_threshold=0.8):
    # TODO: sanity check this.
    # keep the coloc with the highest posterior probability
    logging.info('{:,} pairs tested for coloc'.format(len(coloc)))
    coloced = coloc[coloc['PP.H4.abf']>=PP4_threshold]
    logging.info('{:,} pairs coloc @ PP4 threshold = {}'.format(len(coloced), PP4_threshold))
    #print(coloced.groupby(['xqtl1_modality', 'xqtl1_tissue', 'xqtl1_maf', 'xqtl1_ancestry', 'xqtl1_phenotype', 'xqtl1_cs', 'xqtl2_modality', 'xqtl2_tissue', 'xqtl2_maf', 'xqtl2_ancestry', 'xqtl2_gene']).size().value_counts())
    logging.info('Filtering cases where a single xQTL CS coloced with more than one xQTL CS for the same modality+tissue+maf+ancestry+gene')
    coloced = coloced.sort_values(['PP.H4.abf'], ascending=False).groupby(['xqtl1_modality', 'xqtl1_tissue', 'xqtl1_maf', 'xqtl1_ancestry', 'xqtl1_phenotype', 'xqtl1_cs', 'xqtl2_modality', 'xqtl2_tissue', 'xqtl2_maf', 'xqtl2_ancestry', 'xqtl2_gene']).head(1)
    #print(coloced.groupby(['xqtl1_modality', 'xqtl1_tissue', 'xqtl1_maf', 'xqtl1_ancestry', 'xqtl1_phenotype', 'xqtl1_cs', 'xqtl2_modality', 'xqtl2_tissue', 'xqtl2_maf', 'xqtl2_ancestry', 'xqtl2_gene']).size().value_counts())

    #print(coloced.groupby(['xqtl2_modality', 'xqtl2_tissue', 'xqtl2_maf', 'xqtl2_ancestry', 'xqtl2_phenotype', 'xqtl2_cs', 'xqtl1_modality', 'xqtl1_tissue', 'xqtl1_maf', 'xqtl1_ancestry', 'xqtl1_gene']).size().value_counts())
    coloced = coloced.sort_values(['PP.H4.abf'], ascending=False).groupby(['xqtl2_modality', 'xqtl2_tissue', 'xqtl2_maf', 'xqtl2_ancestry', 'xqtl2_phenotype', 'xqtl2_cs', 'xqtl1_modality', 'xqtl1_tissue', 'xqtl1_maf', 'xqtl1_ancestry', 'xqtl1_gene']).head(1)
    #print(coloced.groupby(['xqtl2_modality', 'xqtl2_tissue', 'xqtl2_maf', 'xqtl2_ancestry', 'xqtl2_phenotype', 'xqtl2_cs', 'xqtl1_modality', 'xqtl1_tissue', 'xqtl1_maf', 'xqtl1_ancestry', 'xqtl1_gene']).size().value_counts())

    logging.info('After filtering, {:,} coloced pairs remain'.format(len(coloced)))
    return coloced


if args.colocs_only:
    coloced = handle_multicoloc(coloc)
    coloced.to_csv(sys.stdout, sep='\t', index=False)
else:
    coloc.to_csv(sys.stdout, sep='\t', index=False)

