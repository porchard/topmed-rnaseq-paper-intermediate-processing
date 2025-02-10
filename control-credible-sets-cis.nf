#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// try per-chrom first
MAF_GLOB = params.maf_glob // {tissue}.{chrom}*
LD_GLOB = params.ld_glob // {tissue}.{chrom}.*
PHENOTYPE_GLOB = params.phenotype_glob // {tissue}.*
SUSIE_GLOB = params.susie_glob // {tissue}.*


process get_controls {

    //container 'docker.io/porchard/general:20220406125608'
    container 'docker.io/porchard/general:20241111'
    memory '30 GB'
    publishDir "${params.results}/controls"

    input:
    tuple val(tissue), val(chrom), path(ld), path(maf), path(phenotypes), path(susie)

    output:
    tuple val(tissue), val(chrom), path("${tissue}.${chrom}.*")

    """
    grep -e variant_id -e ${chrom}_ $susie > susie.chrom.txt

    select-control-credible-sets.py \\
        --credible-sets susie.chrom.txt \\
        --maf $maf \\
        --ld $ld \\
        --phenotypes $phenotypes \\
        --prefix ${tissue}.${chrom}.
    """

}



workflow {
    susie = Channel.fromPath(SUSIE_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, susie
    phenotypes = Channel.fromPath(PHENOTYPE_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, phenotypes
    maf = Channel.fromPath(MAF_GLOB).map({it -> [it.getName().tokenize('.')[0], it.getName().tokenize('.')[1], it]}) // tissue, chrom, maf
    ld = Channel.fromPath(LD_GLOB).map({it -> [it.getName().tokenize('.')[0], it.getName().tokenize('.')[1], it]}) // tissue, chrom, ld
    ld.combine(maf, by: [0, 1]).combine(phenotypes, by: 0).combine(susie, by: 0) | get_controls
}