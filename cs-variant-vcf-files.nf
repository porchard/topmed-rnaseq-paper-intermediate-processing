#!/usr/bin/env nextflow

nextflow.enable.dsl=2

REGIONS = params.regions
VCF_GLOB = params.vcf_glob // name: chrom.vcf
VCF_INDEX_GLOB = params.vcf_index_glob // name: chrom.vcf


process subset_vcf {

    container 'library://porchard/default/general:20220107'
    memory '10 GB'
    time '168h'
    publishDir "${params.results}/vcfs-by-chrom"
    cache 'lenient'

    input:
    tuple val(chrom), path("in.vcf.gz"), path("in.vcf.gz.tbi"), path(regions)

    output:
    tuple val(chrom), path("${chrom}.vcf.gz")

    """
    bcftools view --regions-file $regions -Oz -o ${chrom}.vcf.gz in.vcf.gz
    """

}


process index_vcf {

    container 'library://porchard/default/general:20220107'
    time '168h'
    publishDir "${params.results}/vcfs-by-chrom"
    cache 'lenient'

    input:
    tuple val(chrom), path(vcf)

    output:
    path("*.tbi")

    """
    tabix $vcf
    """

}


workflow {
    vcf_in = Channel.fromPath(VCF_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // chrom, vcf
    vcf_index_in = Channel.fromPath(VCF_INDEX_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // chrom, vcf index
    vcf_in.combine(vcf_index_in, by: 0).combine(Channel.fromPath(REGIONS)) | subset_vcf | index_vcf
}
