#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import topmed_manuscript_clean as tm
import re
import pybedtools as bt
import glob
import os
import sys

VARIANT_LIST, SITES_DIR, ROADMAP_DIR, ENSEMBL_REGULATORY_BUILD, PREFIX = sys.argv[1:]

def get_snpeff_ann(x):
    RE = '^ANN=(.*?)$'
    ann = [re.match(RE, i) for i in x.split(';')]
    ann = [i for i in ann if i is not None]
    if len(ann) != 1:
        raise ValueError(f'No match for {x}')
    return ann[0].group(1)


def parse_snpeff_ann(x):
    l = x.split(',')
    fields = ['allele', 'annotation', 'putative_impact', 'gene_name', 'gene_id', 'feature_type', 'feature_id', 'transcript_biotype', 'rank_total', 'hgvs_c', 'hgvs_p', 'cdna_pos_length', 'cds_pos_len', 'protein_pos_len', 'dist_to_feature', 'messages']
    l = [dict(zip(fields, i.split('|'))) for i in l]
    for i in l:
        i['annotation'] = i['annotation'].split('&')
    return l


def variants_and_bed_to_overlap_matrix(variants, bed):
    # variants is a list of variands (need not be sorted)
    # bed is a BED4 pandas dataframe
    overlap_df = pd.DataFrame(0, columns=sorted(bed['name'].unique()), index=variants)
    variants_bed = tm.variants_to_bed(variants)
    overlaps = bt.BedTool().from_dataframe(variants_bed).sort().intersect(bt.BedTool().from_dataframe(bed).sort(), wa=True, wb=True).to_dataframe()
    for variant_id, name in zip(overlaps['name'], overlaps['thickEnd']):
        overlap_df.at[variant_id,name] = 1
    return overlap_df


variants = pd.read_csv(VARIANT_LIST, header=None, names=['variant_id'])
variants['chrom'] = variants.variant_id.str.split('_', expand=True)[0]



# pull annotations
vcf_annotations = pd.concat([tm.get_variant_info_from_vcf(f'{SITES_DIR}/freeze.9b.{chrom}.pass_and_fail.sites.bcf', df.variant_id.to_list()) for chrom, df in variants.groupby('chrom')])

variant_and_annotation = []
for variant, annotations in zip(vcf_annotations.variant_id, vcf_annotations.INFO):
    for d in parse_snpeff_ann(get_snpeff_ann(annotations)):
        for a in d['annotation']:
            variant_and_annotation.append([variant, a])
snpeff = pd.DataFrame(0, index=variants.variant_id.to_list(), columns=list(set([i[1] for i in variant_and_annotation])))
for (variant, annotation) in variant_and_annotation:
    snpeff.at[variant,annotation] = 1
snpeff.to_csv(f'{PREFIX}snpeff.txt', sep='\t')



# overlap with roadmap annotations
roadmap = {}
ROADMAP_FILES = glob.glob(f'{ROADMAP_DIR}/*_15_coreMarks_hg38lift_mnemonics.bed.gz')
for f in ROADMAP_FILES:
    cell_type = os.path.basename(f).split('_')[0]
    roadmap[cell_type] = variants_and_bed_to_overlap_matrix(variants.variant_id.to_list(), pd.read_csv(f, sep='\t', header=None, names=['chrom', 'start', 'end', 'name']))
for (cell_type, df) in roadmap.items():
    df.to_csv(f'{PREFIX}roadmap_{cell_type}.txt', sep='\t')


# overlap with other annotations
regulatory_build = pd.read_csv(ENSEMBL_REGULATORY_BUILD, sep='\t', header=None)
regulatory_build = regulatory_build.iloc[:,[0, 3, 4, 2]]
regulatory_build.columns = ['chrom', 'start', 'end', 'name']
regulatory_build.chrom = 'chr' + regulatory_build.chrom.astype(str)


regulatory_build = variants_and_bed_to_overlap_matrix(variants.variant_id.to_list(), regulatory_build)
regulatory_build.to_csv(f'{PREFIX}regulatory_build.txt', sep='\t')
