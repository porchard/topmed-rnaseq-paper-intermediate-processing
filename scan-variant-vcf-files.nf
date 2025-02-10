#!/usr/bin/env nextflow

nextflow.enable.dsl=2

VCF_GLOB = params.vcf_glob // name: chrom.vcf
VCF_INDEX_GLOB = params.vcf_index_glob // name: chrom.vcf
SAMPLES_GLOB = params.samples_glob // tissue.*

process subset_vcf {

    container 'library://porchard/default/general:20220107'
    memory '10 GB'
    time '168h'
    publishDir "${params.results}/vcfs-by-chrom"
    cache 'lenient'
    tag "${tissue} ${chrom}"

    input:
    tuple val(chrom), path("in.vcf.gz"), path("in.vcf.gz.tbi"), val(tissue), path(samples)

    output:
    path("${tissue}.${chrom}.vcf.gz")

    script:
    maf = tissue == 'Whole_blood' ? '0.001' : '0.01'

    """
    bcftools view --min-af ${maf}:minor --samples-file $samples -Oz -o ${tissue}.${chrom}.vcf.gz in.vcf.gz
    """

}


process index_vcf {

    container 'library://porchard/default/general:20220107'
    time '168h'
    publishDir "${params.results}/vcfs-by-chrom"
    cache 'lenient'

    input:
    path(vcf)

    output:
    path("*.tbi")

    """
    tabix $vcf
    """

}


workflow {
    vcf_in = Channel.fromPath(VCF_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // chrom, vcf
    vcf_index_in = Channel.fromPath(VCF_INDEX_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // chrom, vcf index
    samples = Channel.fromPath(SAMPLES_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, samples
    vcf_in.combine(vcf_index_in, by: 0).combine(samples) | subset_vcf | index_vcf
}
