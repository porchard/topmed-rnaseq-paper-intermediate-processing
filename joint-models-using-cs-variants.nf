#!/usr/bin/env nextflow

nextflow.enable.dsl=2

COVARIATES_GLOB = params.covariates_glob
PHENOTYPES_GLOB = params.phenotypes_glob
CREDIBLE_SETS_GLOB = params.credible_sets_glob
VCF_PATH = params.vcf_path

process run_model {

	publishDir "${params.results}/joint"
	//container 'docker.io/porchard/general:20230411143108'
	container 'docker.io/porchard/general:20241111'
    tag "${tissue} ${maf} ${modality}"
	memory '100 GB'
	time '24h'

	input:
	tuple val(tissue), val(modality), path(phenotypes), path(covariates), val(maf), path(cs)

	output:
	path("${tissue}.${maf}.${modality}.results.tsv")

	"""
	joint-models-using-cs-variants.py --phenotypes $phenotypes --covariates $covariates --vcf-path $VCF_PATH --credible-sets $cs --output ${tissue}.${maf}.${modality}.results.tsv
	"""

}


workflow {
    covariates = Channel.fromPath(COVARIATES_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}) // tissue, modality, covariates
	phenotypes = Channel.fromPath(PHENOTYPES_GLOB).map({it -> it.getName().tokenize('.')[0..1] + [it]}) // tissue, modality, phenotypes
	cs = Channel.fromPath(CREDIBLE_SETS_GLOB).map({it -> [it.getName().tokenize('.')[0], it.getName().tokenize('.')[2], it.getName().tokenize('.')[1], it]}) // tissue, modality, maf, phenotypes
	phenotypes.combine(covariates, by: [0, 1]).combine(cs, by: [0, 1]) | run_model
}
