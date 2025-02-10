#!/usr/bin/env nextflow

nextflow.enable.dsl=2


PLINK_GLOB = params.plink_glob // {chrom}.*
SAMPLE_GLOB = params.sample_glob // {tissue}.*
MIN_MAF = params.min_maf

process calculate_maf {

    container 'docker.io/porchard/general:20220406125608'
    memory '150 GB'
    publishDir "${params.results}/maf-and-ld"

    input:
    tuple val(tissue), path('samples.txt'), val(chrom), path(plink)

    output:
    tuple val(tissue), val(chrom), path("${tissue}.${chrom}.maf.txt")
    tuple val(tissue), val(chrom), path("${tissue}.${chrom}.ld.txt")

    """
    # subset to samples and MAF threshold
    cat samples.txt | perl -pe 's/^/0\\t/' > keep-samples.txt
    plink -bfile $chrom --keep keep-samples.txt --keep-allele-order --make-bed --maf $MIN_MAF --out subset --threads 1
    # calculate MAF
    plink -bfile subset --freq --threads 1
    mv plink.frq ${tissue}.${chrom}.maf.txt
    # calculate LD
    plink -bfile subset --ld-window 9999999 --memory 150000 --ld-window-kb 2000 --threads 1 --ld-window-r2 0.9 --r2
    mv plink.ld ${tissue}.${chrom}.ld.txt
    """

}



workflow {
    samples = Channel.fromPath(SAMPLE_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, samples
    plink = Channel.fromPath(PLINK_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}).groupTuple() // chrom, plink files
    samples.combine(plink) | calculate_maf
}