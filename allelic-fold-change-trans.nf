#!/usr/bin/env nextflow

SAMPLES_GLOB = params.samples_glob
COVARIATES_GLOB = params.covariates_glob
GENE_VARIANT_PAIRS_GLOB = params.gene_variant_pairs_glob
RNA_COUNTS = params.rna_counts
METADATA = params.metadata
PLINK_GLOB = params.plink_glob

nextflow.enable.dsl=2


process afc {

    publishDir "${params.results}/afc"
    memory '50 GB'
    container 'docker.io/porchard/afc:20230824'

    input:
    tuple val(tissue), path(samples), path(covariates), val(maf), path(gene_variant_pairs), path(rna_counts), path(metadata), path(plink)

    output:
    path("${tissue}.${maf}.txt")

    """
    allelic-fold-change.py --tor $samples --rna-counts $rna_counts --metadata $metadata --covariates $covariates --gene-variant-pairs $gene_variant_pairs --plink-beds *.bed --out ${tissue}.${maf}.txt
    """

}


workflow {
    covariates = Channel.fromPath(COVARIATES_GLOB).map({it -> [it.getName().tokenize('.')[0], it]})
    samples = Channel.fromPath(SAMPLES_GLOB).map({it -> [it.getName().tokenize('.')[0], it]})
    gene_variant_pairs = Channel.fromPath(GENE_VARIANT_PAIRS_GLOB).map({it -> [it.getName().tokenize('.')[0], it.getName().tokenize('.')[1], it]})
    metadata = Channel.fromPath(METADATA)
    rna_counts = Channel.fromPath(RNA_COUNTS)
    plink = Channel.fromPath(PLINK_GLOB).toSortedList().map({it -> [it]})

    samples.combine(covariates, by: 0).combine(gene_variant_pairs, by: 0).combine(rna_counts).combine(metadata).combine(plink) | afc
}