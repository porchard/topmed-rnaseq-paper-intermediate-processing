#!/usr/bin/env python
# coding: utf-8

import sys

from cyvcf2 import VCF
import pandas as pd

SITES_VCF = sys.argv[1]


def parse_snpeff_ann(x):
    l = x.split(',')
    fields = ['allele', 'annotation', 'putative_impact', 'gene_name', 'gene_id', 'feature_type', 'feature_id', 'transcript_biotype', 'rank_total', 'hgvs_c', 'hgvs_p', 'cdna_pos_length', 'cds_pos_len', 'protein_pos_len', 'dist_to_feature', 'messages']
    l = [dict(zip(fields, i.split('|'))) for i in l]
    for i in l:
        i['annotation'] = i['annotation'].split('&')
    return l

variant_and_annotation = []

for variant in VCF(SITES_VCF):
    assert(len(variant.ALT) == 1)
    variant_id = '{}_{}_{}_{}'.format(variant.CHROM, variant.POS, variant.REF, variant.ALT[0])
    ann = variant.INFO.get('ANN')
    for d in parse_snpeff_ann(ann):
        for a in d['annotation']:
            variant_and_annotation.append([variant_id, a])

all_variants = list(set([i[0] for i in variant_and_annotation]))
all_annotations = list(set([i[1] for i in variant_and_annotation]))

snpeff = pd.DataFrame(0, index=all_variants, columns=all_annotations)
for (variant, annotation) in variant_and_annotation:
    snpeff.at[variant,annotation] = 1
snpeff.to_csv(sys.stdout, sep='\t', index=True, header=True, index_label='SNP')