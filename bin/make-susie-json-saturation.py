#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import glob
import re
import os
import json
import sys

# xQTL (joint)
xqtl = []
rda_files = glob.glob('work/tensorqtl-to-coloc-in/saturation/*/maf*/results/coloc-in/*/*.rda')
RE = 'work/tensorqtl-to-coloc-in/saturation/(.*)/(maf\d+)/results/coloc-in/.*/(.*).rda'
for rda_file in rda_files:
    modality, maf, tissue_phenotype = re.match(RE, rda_file).groups()
    tissue, phenotype = tissue_phenotype.split('.')[0], '.'.join(tissue_phenotype.split('.')[1:])
    xqtl.append([modality, maf, tissue, phenotype, 'joint', os.path.abspath(rda_file)])
xqtl = pd.DataFrame(xqtl, columns=['modality', 'maf', 'tissue', 'phenotype', 'ancestry', 'susie'])
xqtl.maf = xqtl.maf.map({'maf001': '1%', 'maf0001': '0.1%', 'maf005': '5%'})
xqtl.modality = xqtl.modality.str.replace('-', '')


# get xQTL regions
# not clipping chromosome ends -- add that later
ciseqtl_regions = pd.concat([pd.read_csv(f, sep='\t', usecols=[0, 1, 2, 3]).assign(tissue=os.path.basename(f).split('.')[0]) for f in glob.glob('data/tensorqtl-in/saturation/cis-eqtl/*.phenotypes.bed.gz')])
ciseqtl_regions['region'] = ciseqtl_regions['#chr'] + '_' + (ciseqtl_regions.start-int(1e6)).map(lambda x: max(0, x)).astype(str) + '_' + (ciseqtl_regions.end+int(1e6)).astype(str)
ciseqtl_regions = ciseqtl_regions[['tissue', 'gene_id', 'region']].rename(columns={'gene_id': 'phenotype'}).drop_duplicates()

# cissqtl_regions = pd.concat([pd.read_csv(f, sep='\t', usecols=[0, 1, 2, 3]).assign(tissue=os.path.basename(f).split('.')[0]) for f in glob.glob('data/saturation/tensorqtl-in/cis-sqtl/*.phenotypes.bed.gz')])
# cissqtl_regions['region'] = cissqtl_regions['#chr'] + '_' + (cissqtl_regions.start-int(1e6)).map(lambda x: max(0, x)).astype(str) + '_' + (cissqtl_regions.end+int(1e6)).astype(str)
# cissqtl_regions = cissqtl_regions[['tissue', 'ID', 'region']].rename(columns={'ID': 'phenotype'}).drop_duplicates()


# xqtl_regions = pd.concat([ciseqtl_regions.assign(modality='ciseqtl'), cissqtl_regions.assign(modality='cissqtl')])
xqtl_regions = ciseqtl_regions.assign(modality='ciseqtl')
length_before = len(xqtl)
xqtl = xqtl.merge(xqtl_regions, on=['tissue', 'phenotype', 'modality'])
length_after = len(xqtl)
assert(length_before == length_after)

assert(xqtl.susie.value_counts().max() == 1)

susie = xqtl.set_index('susie').to_dict(orient='index')

json.dump(susie, sys.stdout, sort_keys = True, indent = 4)