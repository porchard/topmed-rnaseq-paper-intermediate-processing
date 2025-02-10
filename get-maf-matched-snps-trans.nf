#!/usr/bin/env nextflow

nextflow.enable.dsl=2


PLINK_GLOB = params.plink_glob // {chrom}.*
TOP_VARIANTS_GLOB = params.top_variants_glob // {tissue}.*
SAMPLE_GLOB = params.sample_glob // {tissue}.*
MIN_MAF = params.min_maf
SNP_MAPPABILITY = params.snp_mappability

process get_variants_meeting_min_maf_threshold_in_samples {

    container 'docker.io/porchard/general:20220406125608'
    memory '8 GB'

    input:
    tuple val(tissue), path('samples.txt'), val(chrom), path(plink)

    output:
    tuple val(tissue), val(chrom), path("${chrom}.genotypes.bed"), path("${chrom}.genotypes.bim"), path("${chrom}.genotypes.fam")

    """
    cat samples.txt | perl -pe 's/^/0\\t/' > keep-samples.txt
    plink -bfile $chrom --output-chr M --keep-allele-order --make-bed --maf $MIN_MAF --out ${chrom}.genotypes --keep keep-samples.txt
    """

}


process remove_snps_with_mappability_below_1 {

    tag "${tissue} ${chrom}"
    container 'docker.io/porchard/general:20220406125608'
    memory '50 GB'

    input:
    tuple val(tissue), val(chrom), path('genotypes.bed'), path('genotypes.bim'), path('genotypes.fam'), path(snp_mappability), path(snp_mappability_index)

    output:
    tuple val(tissue), val(chrom), path("${chrom}.genotypes.bed"), path("${chrom}.genotypes.bim"), path("${chrom}.genotypes.fam")

    """
    cut -f2 genotypes.bim > variants.txt
    get-variant-mappability.py $snp_mappability variants.txt > variant-mappability.txt
    variant-mappability-to-mappable-variants.py variant-mappability.txt > mappable.txt
    plink -bfile genotypes --output-chr M --keep-allele-order --extract mappable.txt --make-bed --out ${chrom}.genotypes
    """

}


process remove_monomorphic {

    tag "${tissue} ${chrom}"
    container 'docker.io/porchard/general:20220406125608'
    memory '50 GB'

    input:
    tuple val(tissue), val(chrom), path('genotypes.bed'), path('genotypes.bim'), path('genotypes.fam')

    output:
    tuple val(tissue), val(chrom), path("${chrom}.genotypes.*")

    """
    list-monomorphic-fast.py genotypes.bed > drop-variants.txt
    plink -bfile genotypes --output-chr M --keep-allele-order --exclude drop-variants.txt --make-bed --out ${chrom}.genotypes
    """

}


process calculate_maf {

    container 'docker.io/porchard/general:20220406125608'
    memory '30 GB'
    publishDir "${params.results}/maf/${tissue}"

    input:
    tuple val(tissue), val(chrom), path(plink)

    output:
    tuple val(tissue), val(chrom), path("${chrom}.maf.txt")

    """
    plink -bfile ${chrom}.genotypes --freq
    mv plink.frq ${chrom}.maf.txt
    """

}

process get_control_variants {

    container 'docker.io/porchard/general:20220406125608'
    memory '10 GB'
    tag "${tissue} ${chrom}"

    input:
    tuple val(tissue), val(chrom), path(maf), path(top_cs_variants)

    output:
    tuple val(tissue), path("${tissue}.${chrom}.txt")

    """
    # grep ${chrom}_ $top_cs_variants > variants-of-interest.txt
    cat $top_cs_variants | awk '\$1 ~ /${chrom}_/' > variants-of-interest.txt
    n=\$(wc -l < variants-of-interest.txt)
    if [ "\$n" -gt 0  ]
        then
            get-maf-matched-variants.py --plink-maf $maf --min-maf $MIN_MAF --max-maf 0.5 --number-maf-bins 50 --number-controls 5 --variants-of-interest variants-of-interest.txt > ${tissue}.${chrom}.txt
        else
            touch ${tissue}.${chrom}.txt
    fi
    """

}


process merge_chroms {

    publishDir "${params.results}/controls"
    tag "${tissue}"

    input:
    tuple val(tissue), path(x)

    output:
    tuple val(tissue), path("${tissue}.controls.txt")

    """
    cat ${x.join(' ')} | awk 'NR==1' > header.txt
    cp header.txt ${tissue}.controls.txt
    cat ${x.join(' ')} | grep -v -f header.txt >> ${tissue}.controls.txt
    rm header.txt
    """

}



workflow {
    samples = Channel.fromPath(SAMPLE_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, samples
    plink = Channel.fromPath(PLINK_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}).groupTuple() // chrom, plink files
    snp_mappability = Channel.fromPath(SNP_MAPPABILITY)
    snp_mappability_index = Channel.fromPath(SNP_MAPPABILITY + '.tbi')

    variants_meeting_min_maf_threshold = get_variants_meeting_min_maf_threshold_in_samples(samples.combine(plink))
    variants_meeting_min_maf_threshold_mappable_only_no_monomorphic = variants_meeting_min_maf_threshold.combine(snp_mappability).combine(snp_mappability_index) | remove_snps_with_mappability_below_1 | remove_monomorphic
    mafs = calculate_maf(variants_meeting_min_maf_threshold_mappable_only_no_monomorphic) // tissue, chrom, maf

    top_pip_snps = Channel.fromPath(TOP_VARIANTS_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, top variants
    get_control_variants(mafs.combine(top_pip_snps, by: 0)).groupTuple() | merge_chroms

}