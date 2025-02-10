#!/usr/bin/env nextflow

nextflow.enable.dsl=2

SCAN_VARIANT_VCF_GLOB = params.scan_variant_vcf_glob // name: tissue.chrom.vcf
SITES_VCF_GLOB = params.sites_vcf_glob
SITES_VCF_INDEX_GLOB = params.sites_vcf_index_glob

process get_variants {

    container 'library://porchard/default/general:20220107'
    memory '5 GB'
    time '10h'
    publishDir "${params.results}/vcfs-by-chrom"
    cache 'lenient'

    input:
    path(vcf)

    output:
    path("variants.txt")

    """
    bcftools view --drop-genotypes --no-header $vcf | cut -f3 > variants.txt
    """

}


process variants_to_regions {
    
    container 'library://porchard/default/general:20220107'
    memory '25 GB'
    cpus 10
    time '24h'
    publishDir "${params.results}/regions"
    cache 'lenient'

    input:
    path("variants/variants_*.txt")

    output:
    path("variants.bed")

    """
    cat variants/* | sort --parallel=10 -S 20G | uniq | perl -pe 's/_/\\t/g' | awk '{print(\$1, \$2-1, \$2)}' | perl -pe 's/ /\\t/g' > variants.bed
    """

}


process subset_vcf {

    container 'library://porchard/default/general:20220107'
    memory '10 GB'
    time '168h'
    publishDir "${params.results}/sites-by-chrom"
    cache 'lenient'

    input:
    tuple val(chrom), path("in.bcf"), path("in.bcf.csi"), path(regions)

    output:
    tuple val(chrom), path("${chrom}.vcf.gz")

    """
    bcftools view --regions-file $regions -Oz -o ${chrom}.vcf.gz in.bcf
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

process make_annotation_matrices {

    container 'docker.io/porchard/general:20241111'
    time '168h'
    publishDir "${params.results}/annotation-matrices"
    cache 'lenient'

    input:
    tuple val(chrom), path(vcf)

    output:
    path("${chrom}.variant-annotation-matrix.tsv")

    """
    sites-vcf-to-variant-annotation-matrix.py $vcf > ${chrom}.variant-annotation-matrix.tsv
    """

}

process concat_annotation_matrices {

    container 'docker.io/porchard/general:20241111'
    time '168h'
    publishDir "${params.results}/annotation-matrices"
    cache 'lenient'

    input:
    path("matrices/*")

    output:
    path("all.variant-annotation-matrix.tsv")

    """
    concat-variant-annotation-matrices.py matrices/* > all.variant-annotation-matrix.tsv
    """

}



workflow {

    regions = (Channel.fromPath(SCAN_VARIANT_VCF_GLOB) | get_variants).toSortedList() | variants_to_regions
    vcf_in = Channel.fromPath(SITES_VCF_GLOB).map({it -> [it.getName().tokenize('.')[2], it]}) // chrom, vcf
    vcf_index_in = Channel.fromPath(SITES_VCF_INDEX_GLOB).map({it -> [it.getName().tokenize('.')[2], it]}) // chrom, vcf index

    subsetted = vcf_in.combine(vcf_index_in, by: 0).combine(regions) | subset_vcf
    index_vcf(subsetted)
    make_annotation_matrices(subsetted).toSortedList() | concat_annotation_matrices

}
