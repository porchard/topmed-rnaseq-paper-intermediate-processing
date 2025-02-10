#!/usr/bin/env nextflow

nextflow.enable.dsl=2

VARIANT_LIST = params.variants
// CHROMS = (1..22).collect({it -> 'chr' + it}) //  + ['chrX']
PLINK_GLOB = params.plink_glob // {chrom}.*
SAMPLE_GLOB = params.sample_glob // {tissue}.*

process calculate_maf {

    container 'docker.io/porchard/general:20220406125608'
    memory '150 GB'
    publishDir "${params.results}/ld"
    // executor 'local'

    input:
    tuple val(tissue), path('samples.txt'), val(chrom), path(plink), path(variants)

    output:
    tuple val(tissue), val(chrom), path("${tissue}.${chrom}.ld.txt")

    """
    grep ${chrom}_ $variants > keep-variants.${chrom}.txt
    # subset to samples and MAF threshold
    cat samples.txt | perl -pe 's/^/0\\t/' > keep-samples.txt
    plink -bfile $chrom --keep keep-samples.txt --keep-allele-order --make-bed --out subset --threads 1
    # calculate LD
    plink -bfile subset --ld-snp-list keep-variants.${chrom}.txt --ld-window 9999999 --memory 150000 --ld-window-kb 2000 --threads 1 --ld-window-r2 0.8 --r2
    mv plink.ld ${tissue}.${chrom}.ld.txt
    """

}



workflow {
    samples = Channel.fromPath(SAMPLE_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, samples
    plink = Channel.fromPath(PLINK_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}).groupTuple() // chrom, plink files
    samples.combine(plink).combine(Channel.fromPath(VARIANT_LIST)) | calculate_maf
}