#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process replicate {

    container 'docker.io/porchard/general:20241111'
    memory '100 GB'
    publishDir "${params.results}/replicated"

    input:
    path(phenotypes)
    path(covariates)
    val(vcf_dir)
    path(pairs)

    output:
    path('results.txt')

    """
    replicate-trans-in-topmed.py $phenotypes $covariates $vcf_dir $pairs > results.txt
    """

}


workflow {
    replicate(Channel.fromPath(params.phenotypes), Channel.fromPath(params.covariates), params.vcf_dir, Channel.fromPath(params.pairs))
}