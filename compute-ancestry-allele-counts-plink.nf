#!/usr/bin/env nextflow

nextflow.enable.dsl=2


PLINK_GLOB = params.plink_glob // {chrom}.*
NOMINAL_GLOB = params.nominal_glob // {tissue}.*
SAMPLES_GLOB = params.samples_glob // {tissue}.ancestry.*

process get_tested_variants {

    container 'docker.io/porchard/general:20220406125608'
    memory '8 GB'
    publishDir "${params.results}/tested-variants"

    input:
    tuple val(tissue), val(maf), path(nominal)

    output:
    tuple val(tissue), val(maf), path("${tissue}.${maf}.variants.txt")

    """
    zcat $nominal | awk 'NR>1' | cut -f5 | uniq | sort | uniq > "${tissue}.${maf}.variants.txt"
    """

}


process calculate_mac {

    container 'docker.io/porchard/general:20220406125608'
    memory '30 GB'
    publishDir "${params.results}/mac-per-chrom"

    input:
    tuple val(tissue), val(maf), path(variants), val(ancestry), path('samples.txt'), val(chrom), path(plink)

    output:
    tuple val(tissue), val(ancestry), val(maf), val(chrom), path("${tissue}.${ancestry}.${maf}.${chrom}.mac.txt")

    """
    cut -f2 samples.txt | perl -pe 's/^/0\\t/' > keep-samples.txt
    plink -bfile $chrom --freq counts --keep keep-samples.txt --extract $variants
    mv plink.frq.counts ${tissue}.${ancestry}.${maf}.${chrom}.mac.txt
    """

}


process merge_chroms {

    publishDir "${params.results}/mac-merged"
    memory '20 GB'

    input:
    tuple val(tissue), val(ancestry), val(maf), val(chrom), path(x)

    output:
    path("${tissue}.${ancestry}.${maf}.mac.txt")

    """
    cat ${x[0]} | awk 'NR==1' > header.txt
    cp header.txt ${tissue}.${ancestry}.${maf}.mac.txt
    cat ${x.join(' ')} | grep -v -f header.txt >> ${tissue}.${ancestry}.${maf}.mac.txt
    rm header.txt
    """

}



workflow {
    samples = Channel.fromPath(SAMPLES_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}) // tissue, ancestry, samples
    nominal = Channel.fromPath(NOMINAL_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}) // tissue, maf, nominal file
    plink = Channel.fromPath(PLINK_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}).groupTuple() // chrom, plink files
    
    tested_variants = get_tested_variants(nominal) // tissue, maf, variants
    (tested_variants.combine(samples, by: 0).combine(plink) | calculate_mac).groupTuple(by: [0, 1, 2]) | merge_chroms
}