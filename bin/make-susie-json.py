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
rda_files = glob.glob('work/tensorqtl-to-coloc-in/joint/*/maf*/results/coloc-in/*/*.rda')
RE = 'work/tensorqtl-to-coloc-in/joint/(.*)/(maf\d+)/results/coloc-in/.*/(.*).rda'
for rda_file in rda_files:
    modality, maf, tissue_phenotype = re.match(RE, rda_file).groups()
    tissue, phenotype = tissue_phenotype.split('.')[0], '.'.join(tissue_phenotype.split('.')[1:])
    xqtl.append([modality, maf, tissue, phenotype, 'joint', os.path.abspath(rda_file)])
xqtl = pd.DataFrame(xqtl, columns=['modality', 'maf', 'tissue', 'phenotype', 'ancestry', 'susie'])
xqtl.maf = xqtl.maf.map({'maf001': '1%', 'maf0001': '0.1%', 'maf005': '5%'})
xqtl.modality = xqtl.modality.str.replace('-', '')



# get xQTL regions
# not clipping chromosome ends -- add that later
ciseqtl_regions = pd.concat([pd.read_csv(f, sep='\t', usecols=[0, 1, 2, 3]).assign(tissue=os.path.basename(f).split('.')[0]) for f in glob.glob('data/freeze-1.1RNA/files/cis-eqtl/tensorqtl-input/*.phenotypes.bed.gz')])
ciseqtl_regions['region'] = ciseqtl_regions['#chr'] + '_' + (ciseqtl_regions.start-int(1e6)).map(lambda x: max(0, x)).astype(str) + '_' + (ciseqtl_regions.end+int(1e6)).astype(str)
ciseqtl_regions = ciseqtl_regions[['tissue', 'gene_id', 'region']].rename(columns={'gene_id': 'phenotype'}).drop_duplicates()

cissqtl_regions = pd.concat([pd.read_csv(f, sep='\t', usecols=[0, 1, 2, 3]).assign(tissue=os.path.basename(f).split('.')[0]) for f in glob.glob('data/freeze-1.1RNA/files/cis-sqtl/tensorqtl-input/*.phenotypes.bed.gz')])
cissqtl_regions['region'] = cissqtl_regions['#chr'] + '_' + (cissqtl_regions.start-int(1e6)).map(lambda x: max(0, x)).astype(str) + '_' + (cissqtl_regions.end+int(1e6)).astype(str)
cissqtl_regions = cissqtl_regions[['tissue', 'ID', 'region']].rename(columns={'ID': 'phenotype'}).drop_duplicates()

transeqtl_regions = pd.concat([pd.read_csv(f, sep='\t', usecols=[0, 1, 2, 3]).assign(tissue=os.path.basename(f).split('.')[0]) for f in glob.glob('data/tensorqtl-out/trans-eqtl/maf005/trans-susie-phenotypes/*.bed.gz')])
transeqtl_regions['region'] = transeqtl_regions['#chr'] + '_' + (transeqtl_regions.start-int(1e6)).map(lambda x: max(0, x)).astype(str) + '_' + (transeqtl_regions.end+int(1e6)).astype(str)
transeqtl_regions = transeqtl_regions[['tissue', 'gene_id', 'region']].rename(columns={'gene_id': 'phenotype'}).drop_duplicates()

transsqtl_regions = pd.concat([pd.read_csv(f, sep='\t', usecols=[0, 1, 2, 3]).assign(tissue=os.path.basename(f).split('.')[0]) for f in glob.glob('data/tensorqtl-out/trans-sqtl/maf005/trans-susie-phenotypes/*.bed.gz')])
transsqtl_regions['region'] = transsqtl_regions['#chr'] + '_' + (transsqtl_regions.start-int(1e6)).map(lambda x: max(0, x)).astype(str) + '_' + (transsqtl_regions.end+int(1e6)).astype(str)
transsqtl_regions = transsqtl_regions[['tissue', 'gene_id', 'region']].rename(columns={'gene_id': 'phenotype'}).drop_duplicates()


xqtl_regions = pd.concat([ciseqtl_regions.assign(modality='ciseqtl'), cissqtl_regions.assign(modality='cissqtl'), transeqtl_regions.assign(modality='transeqtl'), transsqtl_regions.assign(modality='transsqtl')])
length_before = len(xqtl)
xqtl = xqtl.merge(xqtl_regions, on=['tissue', 'phenotype', 'modality'])
length_after = len(xqtl)
assert(length_before == length_after)


# xQTL (ancestry-specific)
xqtl_ancestry_specific = []
rda_files = glob.glob('work/tensorqtl-to-coloc-in/ancestry-specific/*/maf*/results/coloc-in/*/*.rda')
RE = 'work/tensorqtl-to-coloc-in/ancestry-specific/(.*)/(maf\d+)/results/coloc-in/.*/(.*).rda'
for rda_file in rda_files:
    modality, maf, tissue_phenotype = re.match(RE, rda_file).groups()
    tissue, phenotype = tissue_phenotype.split('.')[0], '.'.join(tissue_phenotype.split('.')[1:])
    tissue, ancestry = tissue.split('___')
    xqtl_ancestry_specific.append([modality, maf, tissue, phenotype, ancestry, os.path.abspath(rda_file)])
xqtl_ancestry_specific = pd.DataFrame(xqtl_ancestry_specific, columns=['modality', 'maf', 'tissue', 'phenotype', 'ancestry', 'susie'])
xqtl_ancestry_specific.maf = xqtl_ancestry_specific.maf.map({'maf001': '1%', 'maf0001': '0.1%', 'maf005': '5%'})
xqtl_ancestry_specific.modality = xqtl_ancestry_specific.modality.str.replace('-', '')



# get xQTL regions
# TODO: not clipping chromosome ends
ciseqtl_regions = pd.concat([pd.read_csv(f, sep='\t', usecols=[0, 1, 2, 3]).assign(tissue=os.path.basename(f).split('.')[0]) for f in glob.glob('data/ancestry-specific/tensorqtl-in/eqtl/results/tensorqtl-in/*.phenotypes.bed.gz')])
ciseqtl_regions['region'] = ciseqtl_regions['#chr'] + '_' + (ciseqtl_regions.start-int(1e6)).map(lambda x: max(0, x)).astype(str) + '_' + (ciseqtl_regions.end+int(1e6)).astype(str)
ciseqtl_regions = ciseqtl_regions[['tissue', 'gene_id', 'region']].rename(columns={'gene_id': 'phenotype'}).drop_duplicates()

cissqtl_regions = pd.concat([pd.read_csv(f, sep='\t', usecols=[0, 1, 2, 3]).assign(tissue=os.path.basename(f).split('.')[0]) for f in glob.glob('data/ancestry-specific/tensorqtl-in/sqtl/results/tensorqtl-in/*.phenotypes.bed.gz')])
cissqtl_regions['region'] = cissqtl_regions['#chr'] + '_' + (cissqtl_regions.start-int(1e6)).map(lambda x: max(0, x)).astype(str) + '_' + (cissqtl_regions.end+int(1e6)).astype(str)
cissqtl_regions = cissqtl_regions[['tissue', 'ID', 'region']].rename(columns={'ID': 'phenotype'}).drop_duplicates()

transeqtl_regions = pd.concat([pd.read_csv(f, sep='\t', usecols=[0, 1, 2, 3]).assign(tissue=os.path.basename(f).split('.')[0]) for f in glob.glob('data/ancestry-specific/tensorqtl-out/trans-eqtl/maf005/trans-susie-phenotypes/*.bed.gz')])
transeqtl_regions['region'] = transeqtl_regions['#chr'] + '_' + (transeqtl_regions.start-int(1e6)).map(lambda x: max(0, x)).astype(str) + '_' + (transeqtl_regions.end+int(1e6)).astype(str)
transeqtl_regions = transeqtl_regions[['tissue', 'gene_id', 'region']].rename(columns={'gene_id': 'phenotype'}).drop_duplicates()

transsqtl_regions = pd.concat([pd.read_csv(f, sep='\t', usecols=[0, 1, 2, 3]).assign(tissue=os.path.basename(f).split('.')[0]) for f in glob.glob('data/ancestry-specific/tensorqtl-out/trans-sqtl/maf005/trans-susie-phenotypes/*.bed.gz')])
transsqtl_regions['region'] = transsqtl_regions['#chr'] + '_' + (transsqtl_regions.start-int(1e6)).map(lambda x: max(0, x)).astype(str) + '_' + (transsqtl_regions.end+int(1e6)).astype(str)
transsqtl_regions = transsqtl_regions[['tissue', 'gene_id', 'region']].rename(columns={'gene_id': 'phenotype'}).drop_duplicates()

xqtl_ancestry_specific_regions = pd.concat([ciseqtl_regions.assign(modality='ciseqtl'), cissqtl_regions.assign(modality='cissqtl'), transeqtl_regions.assign(modality='transeqtl'), transsqtl_regions.assign(modality='transsqtl')])
xqtl_ancestry_specific_regions[['tissue', 'ancestry']] = xqtl_ancestry_specific_regions.tissue.str.split('___', expand=True)



xqtl_ancestry_specific_regions = pd.concat([ciseqtl_regions.assign(modality='ciseqtl'), cissqtl_regions.assign(modality='cissqtl'), transeqtl_regions.assign(modality='transeqtl'), transsqtl_regions.assign(modality='transsqtl')])
xqtl_ancestry_specific_regions[['tissue', 'ancestry']] = xqtl_ancestry_specific_regions.tissue.str.split('___', expand=True)
length_before = len(xqtl_ancestry_specific)
xqtl_ancestry_specific = xqtl_ancestry_specific.merge(xqtl_ancestry_specific_regions, on=['tissue', 'ancestry', 'phenotype', 'modality'])
length_after = len(xqtl_ancestry_specific)
assert(length_before == length_after)



# GWAS
gwas = []
rda_files = glob.glob('data/panukbb-finemapping/work/ancestry-specific-finemapping/lift-susie/results/susie-hg38/*.rda')
RE = 'data/panukbb-finemapping/work/ancestry-specific-finemapping/lift-susie/results/susie-hg38/(.*)___(.*)___(.*).rda'
for rda_file in rda_files:
    trait, ancestry, region = re.match(RE, rda_file).groups()
    gwas.append(['gwas', trait, ancestry, region, os.path.abspath(rda_file)])
gwas = pd.DataFrame(gwas, columns=['modality', 'trait', 'ancestry', 'region', 'susie'])


assert(xqtl.susie.value_counts().max() == 1)
assert(xqtl_ancestry_specific.susie.value_counts().max() == 1)
assert(gwas.susie.value_counts().max() == 1)



susie = gwas.set_index('susie').to_dict(orient='index')
susie.update(xqtl.set_index('susie').to_dict(orient='index'))
susie.update(xqtl_ancestry_specific.set_index('susie').to_dict(orient='index'))

assert(len(susie) == len(gwas) + len(xqtl) + len(xqtl_ancestry_specific))

json.dump(susie, sys.stdout, sort_keys = True, indent = 4)