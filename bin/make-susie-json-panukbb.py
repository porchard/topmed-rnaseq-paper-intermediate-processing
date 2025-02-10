#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import glob
import re
import os
import json
import sys


# GWAS
gwas = []
GLOB = '/net/topmed11/working/porchard/panukbb-finemapping/work/ancestry-specific-finemapping/lift-susie/results/susie-hg38/*.rda'
rda_files = glob.glob(GLOB)
RE = os.path.join(os.path.dirname(GLOB), '(.*)___(.*)___(.*).rda')
for rda_file in rda_files:
    trait, ancestry, region = re.match(RE, rda_file).groups()
    gwas.append(['gwas', trait, ancestry, region, os.path.abspath(rda_file)])
gwas = pd.DataFrame(gwas, columns=['modality', 'trait', 'ancestry', 'region', 'susie'])

assert(gwas.susie.value_counts().max() == 1)
susie = gwas.set_index('susie').to_dict(orient='index')
json.dump(susie, sys.stdout, sort_keys = True, indent = 4)