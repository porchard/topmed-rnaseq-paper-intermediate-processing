#!/usr/bin/env python
# coding: utf-8

import json
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, help='JSON file to filter')
parser.add_argument('--modality', type=str, choices=['ciseqtl', 'cissqtl', 'transeqtl', 'transsqtl', 'gwas'], nargs='*', default=None, help='Modalities to keep')
parser.add_argument('--maf', type=str, choices=['1%', '0.1%', '5%'], nargs='*', default=None, help='MAFs to keep')
parser.add_argument('--tissue', type=str, nargs='*', default=None, help='Tissues to keep')
parser.add_argument('--ancestry', type=str, nargs='*', default=None, choices=['EUR', 'AFR', 'EAS', 'joint'], help='Ancestries to keep')
parser.add_argument('--trait', type=str, nargs='*', default=None, help='Traits to keep (GWAS only)')
parser.add_argument('--phenotype', type=str, nargs='*', default=None, help='Phenotypes (genes/splicing phenotypes) to keep (xQTL only)')
args = parser.parse_args()

with open(args.input) as f:
    data = json.load(f)

filtered = dict()

for susie, info in data.items():
    if args.modality is not None:
        if 'modality' in info and info['modality'] not in args.modality:
            continue
    if args.maf is not None:
        if 'maf' in info and info['maf'] not in args.maf:
            continue
    if args.tissue is not None:
        if 'tissue' in info and info['tissue'] not in args.tissue:
            continue
    if args.ancestry is not None:
        if 'ancestry' in info and info['ancestry'] not in args.ancestry:
            continue
    if args.trait is not None:
        if 'trait' in info and info['trait'] not in args.trait:
            continue
    if args.phenotype is not None:
        if 'phenotype' in info and info['phenotype'] not in args.phenotype:
            continue
    filtered[susie] = info

logging.info(f'Filtered {len(data)} to {len(filtered)}')

print(json.dumps(filtered, indent=4, sort_keys=True))
