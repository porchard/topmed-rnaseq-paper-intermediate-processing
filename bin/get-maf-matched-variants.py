#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import logging
import sys
import argparse
import glob
import os
import logging
import argparse

parser = argparse.ArgumentParser('')
parser.add_argument('--min-maf', type=float, help='Minimum MAF to use for setting MAF bins')
parser.add_argument('--max-maf', type=float, help='Maximum MAF to use for setting MAF bins')
parser.add_argument('--number-maf-bins', type=int, help='Number of MAF bins to use')
parser.add_argument('--number-controls', type=int, help='Number of control variants to select.')
parser.add_argument('--plink-maf', help='Plink MAF output', nargs='+')
parser.add_argument('--variants-of-interest', help='File listing variants to select controls for. Variants may be listed multiple times, in which case N * number_controls control variants will be fetched for it.')
args = parser.parse_args()


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

# options:
# number of variants
# number of MAF bins
# assume that any filtering (on samples, on variants for mappability etc) has been done upstream
# in that case, just need a table of SNP --> MAF, and a list of the SNPs we're trying to find matches for
# NUMBER_MAF_BINS = 3
# MIN_MAF = 0.01
# MAX_MAF = 0.50
# NUMBER_CONTROLS = 10
# SELECT_CONTROLS_FOR = ['chr22_15306199_A_T', 'chr22_15319953_C_T']

NUMBER_MAF_BINS = args.number_maf_bins
MIN_MAF = args.min_maf
MAX_MAF = args.max_maf
NUMBER_CONTROLS = args.number_controls
SELECT_CONTROLS_FOR = []
with open(args.variants_of_interest, 'r') as fh:
    for line in fh:
        SELECT_CONTROLS_FOR.append(line.rstrip())



STEP = (MAX_MAF-MIN_MAF) / NUMBER_MAF_BINS
MAF_BINS = list(np.arange(start=MIN_MAF, stop=MAX_MAF, step=STEP)) + [MAX_MAF + 0.0001]
logging.info('Using {} MAF BINS. Boundaries: {}'.format(len(MAF_BINS) - 1, ', '.join([str(i) for i in MAF_BINS])))


mafs = pd.concat([pd.read_csv(f, delim_whitespace=True) for f in args.plink_maf])
mafs['maf_bin'] = np.digitize(mafs.MAF, MAF_BINS)
variant_to_maf_bin = dict(zip(mafs.SNP, mafs.maf_bin))
variant_to_maf = dict(zip(mafs.SNP, mafs.MAF))


if not all(pd.Series(SELECT_CONTROLS_FOR).isin(mafs.SNP)):
    raise ValueError('All of the variants for which youd like controls must be included in the MAF table')


snps_of_interest = mafs.set_index('SNP').loc[SELECT_CONTROLS_FOR].reset_index()
possible_controls = mafs[~mafs.SNP.isin(SELECT_CONTROLS_FOR)]
controls = []


for maf_bin, df in snps_of_interest.groupby('maf_bin'):
    maf_bin_possible_controls = possible_controls[possible_controls.maf_bin==maf_bin]
    if not len(df) * NUMBER_CONTROLS <= len(maf_bin_possible_controls):
        raise ValueError('Not enough potential control SNPs')
    control_variants = maf_bin_possible_controls.sample(n=len(df) * NUMBER_CONTROLS, replace=False).SNP.to_list()
    for count, variant in enumerate(df.SNP.to_list()):
        for control_variant in control_variants[(count*NUMBER_CONTROLS):((1+count)*NUMBER_CONTROLS)]:
            controls.append([variant, control_variant])
controls = pd.DataFrame(controls, columns=['variant', 'control_variant'])
assert(all(controls.variant.map(variant_to_maf_bin) == controls.control_variant.map(variant_to_maf_bin)))
controls['variant_maf'] = controls.variant.map(variant_to_maf)
controls['control_variant_maf'] = controls.control_variant.map(variant_to_maf)


controls.to_csv(sys.stdout, sep='\t', index=False)