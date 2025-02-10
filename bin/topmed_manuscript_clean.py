#!/usr/bin/env python
# coding: utf-8

import glob
import os
import re
import gzip
import tempfile
import subprocess
import sys

import pandas as pd
import numpy as np
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import matplotlib as mpl
import pybedtools as bt
import qtl.norm as norm
import statsmodels.api as sm
from cyvcf2 import VCF


ANCESTRY_ABBREVIATIONS = {
    'EUR': 'Europe',
    'AFR': 'Sub-saharan Africa',
    'AMR': 'Native America',
    'EAS': 'East Asia',
    'MES': 'Middle East',
    'SAS': 'Central and South Asia'
}

cis_eqtl_pcs = {
    'Lung': 75,
    'Monocyte': 30,
    'Nasal_epithelial': 30,
    'PBMC': 30,
    'T_cell': 30,
    'Whole_blood': 100,
}

cis_sqtl_pcs = {
    'Lung': 10,
    'Monocyte': 10,
    'Nasal_epithelial': 10,
    'PBMC': 10,
    'T_cell': 10,
    'Whole_blood': 10,
}

topmed_to_roadmap = {
    'Whole_blood': 'E062',
    'PBMC': 'E062',
    'Lung': 'E096',
    'Monocyte': 'E029',
    'T_cell': 'E034',
    'Nasal_epithelial': 'E114'
}

roadmap_state_names = {
    '1_TssA': 'Active_TSS',
    '2_TssAFlnk': 'Flanking Active TSS',
    '3_TxFlnk':	"Transcr. at gene 5' and 3'",
    '4_Tx':	'Strong transcription',
    '5_TxWk': 'Weak transcription',
    '6_EnhG': 'Genic enhancers',
    '7_Enh': 'Enhancers',
    '8_ZNF/Rpts': 'ZNF genes & repeats',
    '9_Het': 'Heterochromatin',
    '10_TssBiv': 'Bivalent/Poised TSS',
    '11_BivFlnk': 'Flanking Bivalent TSS/Enh',
    '12_EnhBiv': 'Bivalent Enhancer',
    '13_ReprPC': 'Repressed PolyComb',
    '14_ReprPCWk': 'Weak Repressed PolyComb',
    '15_Quies': 'Quiescent/Low'
}

tissues = ['Lung', 'Monocyte', 'Nasal_epithelial', 'PBMC', 'T_cell', 'Whole_blood']
modality = ['cis-eQTL', 'cis-sQTL', 'trans-eQTL', 'trans-sQTL']

palettes = {
    'tissue' : {
        'Whole_blood': '#1b9e77',
        'Nasal_epithelial': '#d95f02',
        'T_cell': '#7570b3',
        'Monocyte': '#e7298a',
        'PBMC': '#66a61e',
        'Lung': '#e6ab02'
    },
    'modality': {'cis-eQTL': '#0343df', 'cis-sQTL': '#15b01a', 'trans-eQTL': '#95d0fc', 'trans-sQTL': '#96f97b'}
}


def make_colormap_dict(keys, palette='viridis'):
    """
    Given list of items (in order), create a dict of item --> color.

    Input:
    keys: list of items.
    palette: name of matplotlib color palette to use

    Returns: Dict of item --> color (hex)
    """
    assert(isinstance(keys, list))
    assert(isinstance(palette, str))
    cmap = cm.get_cmap(palette, len(keys))
    return {keys[i]: cmap(i) for i in range(cmap.N)}


# https://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def make_locuszoom_cmap_and_norm(default=True, include_negative=False, palette=None):
    if default:
        lz_colors = ["#7F7F7F", "#282973", "#8CCCF0", "#69BD45", "#F9A41A", "#ED1F24"]
        bounds = np.append(-1, np.arange(0,1.2,0.2))
        if include_negative:
            cmap = colors.ListedColormap(lz_colors)
            norm = colors.BoundaryNorm(bounds, cmap.N)
        else:
            cmap = colors.ListedColormap(lz_colors[1:])
            norm = colors.BoundaryNorm(bounds[1:], cmap.N)
    else:
        if palette is None:
            raise ValueError('If default = False, must provide palette')
        cmap = cm.get_cmap(palette)
        cmap = truncate_colormap(cmap, 0, 0.7, n=cmap.N)
        bounds = [0, 0.2, 0.4, 0.6, 0.8, 1.0] if not include_negative else [-1, 0, 0.2, 0.4, 0.6, 0.8, 1.0]
        norm = colors.BoundaryNorm(bounds, cmap.N)
    return (cmap, norm)


def make_locuszoom_legend(cmap, norm):
    fig, ax = plt.subplots(figsize=(0.2,2))

    fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                cax=ax, orientation='vertical', label='$r^2$')
    
    return (fig, ax)


def phenotype_id_to_gene_id(x):
    # first, try to match with a version
    ENSEMBL_RE_WITH_VERSION = r'ENSG\d+\.\d+'
    ENSEMBL_RE_WITHOUT_VERSION = r'ENSG\d+'
    with_version = re.search(ENSEMBL_RE_WITH_VERSION, x)
    without_version = re.search(ENSEMBL_RE_WITHOUT_VERSION, x)
    if with_version:
        return with_version.group(0)
    elif without_version:
        return without_version.group(0)
    else:
        raise ValueError(f'Not able to infer gene ID from {x}')


def format_modality(modality):
    modality_lower = modality.lower()
    assert('trans' in modality_lower or 'cis' in modality_lower)
    cis_trans = 'trans' if 'trans' in modality_lower else 'cis'
    assert('eqtl' in modality_lower or 'sqtl' in modality_lower)
    e_sqtl = 'eQTL' if 'eqtl' in modality_lower else 'sQTL'
    return cis_trans + '-' + e_sqtl


def residualize(x, covariates):
    assert(isinstance(x, pd.Series))
    assert(isinstance(covariates, pd.DataFrame))
    model = sm.OLS(x, sm.add_constant(covariates.T)).fit()
    return model.resid


def associate(phenotypes, genotypes, covariates):
    assert(isinstance(phenotypes, pd.Series))
    assert(isinstance(genotypes, pd.Series))
    assert(isinstance(covariates, pd.DataFrame))
    # residualize
    ADD_CONSTANT = True
    p_residuals = residualize(phenotypes, covariates)
    g_residuals = residualize(genotypes, covariates)
    p_residuals = norm.inverse_normal_transform(p_residuals)
    model = sm.OLS(p_residuals, sm.add_constant(g_residuals)) if ADD_CONSTANT else sm.OLS(p_residuals, g_residuals)
    number_samples = len(phenotypes)
    number_covariates = len(covariates)
    dof = number_samples - number_covariates - 2
    model.df_model = dof
    model.df_resid = dof
    results = model.fit()
    return results


def parse_attribute(attribute_series: pd.Series, attribute_name: str) -> pd.Series:
    """
    Parse the attributes column of a (GENCODE/RefSeq) GTF file.

    Input:
    * a [str]: the attributes element (column 9 of the GTF file)
    * regex [str]: a regular expression that will be iteratively applied to the attribute string to capture attribute key, val pairs. Default should work for GENCODE/RefSeq
    """
    if not isinstance(attribute_series, pd.Series):
        raise TypeError('attribute_series must be a pandas Series')
    if not isinstance(attribute_name, str):
        raise TypeError('attribute_name must be a string')
    
    return attribute_series.str.extract(f'{attribute_name} "(.*?)"')


def gtf_to_df(gtf: str, parse_attributes: list=None) -> pd.DataFrame:
    df = pd.read_csv(gtf, sep='\t', header=None, names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'], comment='#')
    if parse_attributes is not None:
        for a in parse_attributes:
            df[a] = parse_attribute(df.attributes, a)
    return df


def get_alt_allele_counts(vcf_path, region=None, variants=None, samples=None):
    if region is not None and variants is not None:
        sys.stderr.write('Variants and regions is set; only fetching variants from that region\n')
    
    vcf = VCF(vcf_path, gts012=True, samples=samples)
    alt_allele_counts = []
    variants_pulled = []

    if region is None and variants is not None:
        for variant in variants:
            chrom, pos = variant.split('_')[:2]
            start, end = int(pos) - 1, int(pos)
            for v in vcf(f'{chrom}:{start}-{end}'):
                if v.ID == variant:
                    alt_allele_counts.append(list(v.gt_types))
                    variants_pulled.append(v.ID)
    elif region is not None:
        if variants is not None:
            variant_set = set(variants)
        for v in vcf(region):
            if variants is not None and v.ID not in variant_set:
                continue
            alt_allele_counts.append(list(v.gt_types))
            variants_pulled.append(v.ID)
    else:
        for v in vcf():
            alt_allele_counts.append(list(v.gt_types))
            variants_pulled.append(v.ID)

    alt_allele_counts_df = pd.DataFrame(alt_allele_counts, columns=vcf.samples, index=variants_pulled)
    return alt_allele_counts_df


def get_genotypes_from_vcf(vcf_path, region=None, variants=None, samples=None):
    """TODO: Handle phasing"""
    df = get_alt_allele_counts(vcf_path, region=region, variants=variants, samples=samples)
    assert(df.isin([0, 1, 2]).all().all())
    return df.replace({0: '0/0', 1: '0/1', 2: '1/1'})


def alt_allele_counts_matrix_to_r2_matrix(x):
    """
    Rows are variants, in format {chrom}_{pos}_{ref}_{alt}
    Columns are samples
    Variants on separate chromosomes are made to have LD = 0
    """

    # split by chrom; all between-chrom values are set to 0
    r2_by_chrom = []
    for chrom, df in x.groupby(lambda x: x.split('_')[0]):
        r2_by_chrom.append(df.T.corr() ** 2)
    return pd.concat(r2_by_chrom).fillna(0)


def alt_allele_counts_matrix_to_r2_series(x, variant):
    """
    Rows are variants, in format {chrom}_{pos}_{ref}_{alt}
    Columns are samples
    Variants on separate chromosomes are made to have LD = 0
    """
    ld_with_variant = x.T.corrwith(x.loc[variant]) ** 2
    mask = (ld_with_variant.index.to_series().str.split('_', expand=True)[0] == variant.split('_')[0])
    ld_with_variant[~mask] = 0
    return ld_with_variant


def parse_region(s):
    chrom, start, end = re.match('^(.*)[-:_](\d+)[-:_](\d+)$', s).groups()
    start, end = int(start), int(end)
    return (chrom, start, end)


def read_vcf_sample_line(vcf_path):
    x = VCF(vcf_path)
    return ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + x.samples


def read_tabix_header(f):
    is_vcf = re.match('^.*\\.vcf(.gz)?$', os.path.basename(f))
    is_bcf = re.match('^.*\\.bcf(.gz)?$', os.path.basename(f))
    if is_vcf or is_bcf:
        header = read_vcf_sample_line(f)
        return header
    else:
        with gzip.open(f, 'rt') as fh:
            header = fh.readline().rstrip().split()
        # if there's no #, it's not actually a header
        if header[0].startswith('#'):
            header[0] = header[0].replace('#', '')
            return header
        else:
            return None


def tabix(f, region):
    """
    f: file to tabix
    region: region, or list of regions, in format chr:start-end (can use :,-,_ as separators)
    """
    if not isinstance(region, list) and not isinstance(region, str):
        raise TypeError('region must be a list or a str')
    if not isinstance(f, str):
        raise TypeError('f must be a str (path to a file)')
    if isinstance(region, str):
        region = [region]
    
    regions = []
    for r in region:
        chrom, start, end = parse_region(r)
        start, end = str(start), str(end)
        regions.append([chrom, start, end])
    with tempfile.NamedTemporaryFile() as tmpf:
        pd.DataFrame(regions).to_csv(tmpf.name, sep='\t', index=False, header=False)
        sp = subprocess.run(['tabix', '--regions', tmpf.name, f], capture_output=True, check=True)
    txt = sp.stdout.decode().split('\n')
    if txt[-1] == '':
        txt = txt[:-1]
    
    if len(txt) == 0:
        return None
    return pd.DataFrame([i.split('\t') for i in txt], columns=read_tabix_header(f))


def variants_to_bed (variants: list) -> pd.DataFrame:
    # variants is a list of variants in form chrom_pos_ref_alt
    # returns dataframe of chrom, start, end, variant
    if not isinstance(variants, (list, pd.Series, np.ndarray)):
        raise ValueError('variants must be a list, series, or array')
    df = pd.Series(variants).str.split('[_:]', expand=True, regex=True) if isinstance(variants, (np.ndarray, list)) else variants.str.split('[_:]', expand=True, regex=True)
    assert(len(df.columns) in [2, 4])
    df.columns = ['chrom', 'pos', 'ref', 'alt'] if len(df.columns) == 4 else ['chrom', 'pos']
    df['end'] = df.pos.astype(int)
    df['start'] = df.end - 1
    df['name'] = list(variants)
    return df[['chrom', 'start', 'end', 'name']]


def bed_to_regions(bed_df: pd.DataFrame, merge: bool=False) -> list:
    if merge:
        tmp = bt.BedTool().from_dataframe(bed_df).sort().merge().to_dataframe().iloc[:,[0, 1, 2]]
        tmp.columns = ['chrom', 'start', 'end']
    else:
        tmp = bed_df[['chrom', 'start', 'end']]
    return list(tmp.chrom.astype(str) + ':' + tmp.start.astype(str) + '-' + tmp.end.astype(str))


def variants_to_regions(x: list):
    # variants can be of format:
    # chrom_pos_.*
    # chrom:pos:.*
    if len(x) == 0:
        return []
    else:
        return bed_to_regions(variants_to_bed(x))


def get_variant_info_from_vcf(vcf_path, variants):
    vcf = VCF(vcf_path, gts012=True, samples=None)

    variants_pulled = []
    ids = []
    info = []

    for variant in variants:
        chrom, pos = variant.split('_')[:2]
        start, end = int(pos) - 1, int(pos)
        for v in vcf(f'{chrom}:{start}-{end}'):
            assert(len(v.ALT) == 1)
            variant_id = v.CHROM + '_' + str(v.POS) + '_' + v.REF + '_' + v.ALT[0]
            if variant_id == variant:
                variants_pulled.append(variant_id)
                ids.append(v.ID)
                info_string = ';'.join(['{}={}'.format(i[0], i[1]) for i in v.INFO])
                info.append(info_string)
    return pd.DataFrame({'variant_id': variants_pulled, 'ID': ids, 'INFO': info})


def ld(v1, v2):
    # note: assumes bi-allelic
    # returns np.nan
    # v1 and v2 are vectors of phased genotypes, e.g. [0|1, 0|0, 0|1]
    # A, B are alt alleles
    v1 = pd.Series(v1).str.split('|', expand=True, regex=False)
    v2 = pd.Series(v2).str.split('|', expand=True, regex=False)
    v1_alt_allele_dosages = v1.astype(int).sum(axis=1)
    v2_alt_allele_dosages = v2.astype(int).sum(axis=1)
    haplotypes = (v1+v2)
    haplotypes = pd.Series(haplotypes[0].to_list() + haplotypes[1].to_list())
    haplotypes = haplotypes.str.split('', expand=True)[[1, 2]].astype(int)
    haplotypes.columns = ['A', 'B']
    pA = haplotypes.A.mean()
    pB = haplotypes.B.mean()
    pAB = ((haplotypes.A==1) & (haplotypes.B==1)).mean()
    pAb = ((haplotypes.A==1) & (haplotypes.B==0)).mean()
    paB = ((haplotypes.A==0) & (haplotypes.B==1)).mean()
    pab = ((haplotypes.A==0) & (haplotypes.B==0)).mean()
    D = pAB - (pA*pB)
    D_max = max([-1*pA*pB, -1*(1-pA)*(1-pB)]) if D < 0 else min([pA*(1-pB), (1-pA)*pB])
    D_prime = D / D_max
    r2 = (D**2) / (pA * (1-pA) * pB * (1-pB))
    r2_alt_allele_dosages = v1_alt_allele_dosages.corr(v2_alt_allele_dosages)**2
    return (D, D_max, D_prime, r2, r2_alt_allele_dosages)


def top_pip_variants(cs):
    """
    Given a dataframe representing CS (having columns ['phenotype_id', 'variant_id', 'cs_id', 'pip']; can have other columns too),
    return a dataframe containing only the top PIP variant per credible set (credible set defined by phenotype_id + csID pair)
    """
    # Validate input
    if not isinstance(cs, pd.DataFrame):
        raise TypeError('cs must be a DataFrame')
    for i in ['phenotype_id', 'variant_id', 'cs_id', 'pip']:
        if not i in cs.columns.to_list():
            raise ValueError(f'cs must include column {i}')

    return cs.sort_values(['pip', 'variant_id'], ascending=[False, True]).groupby(['phenotype_id', 'cs_id']).head(1)


def bin_integers(x, bins):
    """
    x: [0, 5, 3, 1]
    bins = ['0', '1', '2-4', '5+']
    return: ['0', '5+', '2-4', '1']
    """
    mappings = dict()
    for i in bins:
        if '+' in i:
            for j in range(int(i.replace('+', '')), max(x)+1):
                assert(j) not in mappings
                mappings[j] = i
        elif '-' in i:
            for j in range(int(i.split('-')[0]), int(i.split('-')[1])+1):
                assert(j) not in mappings
                mappings[j] = i
        else:
            assert(int(i)) not in mappings
            mappings[int(i)] = i
    y = [mappings[i] for i in x]
    return pd.Categorical(y, categories=bins, ordered=True)

# given tissue and modality, load scan files (phenotypes, covariates)
def load_scan_files(tissue, modality, ancestry='joint'):
    tissues = ['Lung', 'PBMC', 'Monocyte', 'Nasal_epithelial', 'T_cell', 'Whole_blood']
    modalities = ['cis-eQTL', 'cis-sQTL', 'trans-eQTL', 'trans-sQTL']
    assert(modality in modalities)
    assert(tissue in tissues)
    assert(ancestry == 'joint') # not supporting anything else for now

    input_files = glob.glob('/net/topmed10/working/porchard/rnaseq/work/freezes/freeze-1.1RNA/files/cis-*qtl/tensorqtl-input/*')
    INPUT_FILE_RE = re.compile('/net/topmed10/working/porchard/rnaseq/work/freezes/freeze-1.1RNA/files/(cis-.qtl)/tensorqtl-input/(.*).tensorqtl-in.[\d\.]*(phenotypes|covariates)..*')
    scan_files = {}
    
    for f in input_files:
        cis_modality, file_tissue, phenotypes_or_covariates = INPUT_FILE_RE.match(f).groups()
        cis_modality = format_modality(cis_modality)
        if cis_modality == 'cis-eQTL':
            scan_files[('cis-eQTL', file_tissue, phenotypes_or_covariates)] = f
            scan_files[('trans-eQTL', file_tissue, phenotypes_or_covariates)] = f
        elif cis_modality == 'cis-sQTL':
            scan_files[('cis-sQTL', file_tissue, phenotypes_or_covariates)] = f
            scan_files[('trans-sQTL', file_tissue, phenotypes_or_covariates)] = f
    
    phenotypes = pd.read_csv(scan_files[(modality, tissue, 'phenotypes')], sep='\t', index_col=3).drop(columns=['#chr', 'start', 'end'])
    covariates = pd.read_csv(scan_files[(modality, tissue, 'covariates')], sep='\t', index_col=0)
    assert(covariates.columns.to_list() == phenotypes.columns.to_list())

    if modality == 'trans-eQTL' and tissue == 'Whole_blood':
        # remove dropped PCs
        dropped_pcs = [25] + list(range(51, 101))
        covariates = covariates.drop([f'phenotype_PC{i}' for i in dropped_pcs])

    return (phenotypes, covariates)



def compare_credible_sets(scan_1_cs, scan_2_cs, summarize=True):
    """
    Given a dataframe representing scan_1_cs and scan_2_cs (each having columns ['phenotype_id', 'variant_id', 'cs_id']; can have other columns too),
    return a dataframe showing, for each CS in scan_1_cs, whether it overlaps a CS in scan_2_cs (for the same phenotype_id)
    If summarize = False, simply returns scan_1_cs with an added column indicating whether each credible set SNP is a credible set SNP 
    for the same phenotype in the other scan
    """
    # Validate input
    if not isinstance(scan_1_cs, pd.DataFrame):
        raise TypeError('scan_1_cs must be a DataFrame')
    if not isinstance(scan_2_cs, pd.DataFrame):
        raise TypeError('scan_2_cs must be a DataFrame')
    for i in ['phenotype_id', 'variant_id', 'cs_id']:
        if not i in scan_1_cs.columns.to_list():
            raise ValueError(f'scan_1_cs must include column {i}')
        if not i in scan_2_cs.columns.to_list():
            raise ValueError(f'scan_2_cs must include column {i}')

    results = scan_1_cs.merge(scan_2_cs[['phenotype_id', 'variant_id']].drop_duplicates().assign(in_other_scan_cs=1), how='left')
    results.in_other_scan_cs = results.in_other_scan_cs.fillna(0).astype(int)
    assert(len(results) == len(scan_1_cs))

    if summarize:
        return results.groupby(['phenotype_id', 'cs_id']).in_other_scan_cs.max().reset_index()
    else:
        return results



def add_legend_from_colors(d, ax, loc='best', marker='o'):
    """
    Given a dictionary from label --> color and a matplotlib ax,
    add a legend to the ax
    """
    legend_elements = [Line2D([0], [0], marker=marker, color=color, label=label) for label, color in d.items()]
    ax.legend(handles=legend_elements, loc=loc)
    return ax


def rotate_xticklabels(ax, rot=45, ha='right'):
    for t in ax.get_xticklabels():
        t.set(rotation=rot, ha=ha)
    return True


@ticker.FuncFormatter
def bmk_formatter(x, pos):
    """
    Tick label formatting function that converts labels to B/M/k (billions, millions, thousands).

    Usage:
    ax.xaxis.set_major_formatter(xformatter)
    """
    if abs(x) >= 1e9:
        return '{}B'.format(x/1e9)
    elif abs(x) >= 1e6:
        return '{}M'.format(x/1e6)
    elif abs(x) >= 1e3:
        if x % 1e3 == 0:
            return '{}k'.format(round(x/1e3))
        else:
            return '{}k'.format(x/1e3)
    else:
        return x


@ticker.FuncFormatter
def pos_formatter(x, pos):
    """
    Tick label formatting function that converts 1e9 --> Gb, 1e6 --> Mb, 1e3 --> kb.

    Usage:
    ax.xaxis.set_major_formatter(pos_formatter)
    """
    if abs(x) >= 1e9:
        return '{} Gb'.format(x/1e9)
    elif abs(x) >= 1e6:
        return '{} Mb'.format(x/1e6)
    elif abs(x) >= 1e3:
        return '{} kb'.format(x/1e3)
    else:
        return x
