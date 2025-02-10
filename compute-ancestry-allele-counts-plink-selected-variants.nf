#!/usr/bin/env nextflow

nextflow.enable.dsl=2


PLINK_GLOB = params.plink_glob // {chrom}.*
VARIANTS = params.variants
SAMPLES_GLOB = params.samples_glob // {tissue}.ancestry.*

process calculate_mac {

    container 'docker.io/porchard/general:20220406125608'
    memory '30 GB'
    publishDir "${params.results}/mac-per-chrom"

    input:
    tuple val(tissue), val(ancestry), path('samples.txt'), path(variants), val(chrom), path(plink)

    output:
    tuple val(tissue), val(ancestry), val(chrom), path("${tissue}.${ancestry}.${chrom}.mac.txt")

    """
    cut -f2 samples.txt | perl -pe 's/^/0\\t/' > keep-samples.txt
    plink -bfile $chrom --freq counts --keep keep-samples.txt --extract $variants
    mv plink.frq.counts ${tissue}.${ancestry}.${chrom}.mac.txt
    """

}


process merge_chroms {

    publishDir "${params.results}/mac-merged"

    input:
    tuple val(tissue), val(ancestry), val(chrom), path(x)

    output:
    path("${tissue}.${ancestry}.mac.txt")

    """
    cat ${x[0]} | awk 'NR==1' > header.txt
    cp header.txt ${tissue}.${ancestry}.mac.txt
    cat ${x.join(' ')} | grep -v -f header.txt >> ${tissue}.${ancestry}.mac.txt
    rm header.txt
    """

}



workflow {
    samples = Channel.fromPath(SAMPLES_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}) // tissue, ancestry, samples
    variants = Channel.fromPath(VARIANTS)
    plink = Channel.fromPath(PLINK_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}).groupTuple() // chrom, plink files
    
    (samples.combine(variants).combine(plink) | calculate_mac).groupTuple(by: [0, 1]) | merge_chroms
}