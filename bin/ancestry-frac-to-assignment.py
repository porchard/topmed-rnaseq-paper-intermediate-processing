#!/usr/bin/env python
# coding: utf-8

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('metadata')
parser.add_argument('ancestry')
parser.add_argument('--threshold', type=float, default=0.0, help='Min. fraction from the most represented ancestry to be assigned to that ancestry')
args = parser.parse_args()

metadata = pd.read_csv(args.metadata, sep='\t')
metadata = metadata[metadata.used_for_scan]
ancestry = pd.read_csv(args.ancestry, sep='\t', index_col=0)

sample_ancestries = ancestry.loc[metadata.wgs.unique(),:]
max_frac = sample_ancestries.max(axis=1)
max_ancestry = sample_ancestries.idxmax(axis=1)
wgs_to_ancestry = max_ancestry[max_frac>=args.threshold].to_dict()

metadata['ancestry'] = metadata.wgs.map(wgs_to_ancestry)

for i, df in metadata.groupby(['tissue', 'ancestry']):
    tissue, ancestry = i
    df[['tor', 'wgs']].to_csv(f'{tissue}.{ancestry}.samples.txt', sep='\t', index=False, header=False)