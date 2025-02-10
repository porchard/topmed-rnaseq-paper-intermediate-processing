#!/usr/bin/env nextflow

nextflow.enable.dsl=2
PAIRS = params.pairs
LABEL_1 = params.label_1
LABEL_2 = params.label_2

process run_coloc {

    container 'docker.io/porchard/coloc:20221220'
    maxForks 50
    memory '10 GB'

    input:
    path(x)

    output:
    path('coloc.txt')

    """
    run-coloc-batch-preload.R --label-1 $LABEL_1 --label-2 $LABEL_2 $x coloc.txt
    """

}


process concat {

    publishDir "${params.results}/txt"
    beforeScript 'ulimit -Ss unlimited'


    input:
    path("coloc.*.txt")

    output:
    path("coloc.txt")

    """
    cat coloc.1.txt | awk 'NR==1' > header.txt
    cp header.txt coloc.txt
    cat coloc.*.txt | grep -v -f header.txt >> coloc.txt
    rm header.txt
    """
}

workflow {
    coloc = Channel.fromPath(PAIRS).splitText(by: 60000, file: true) | run_coloc
    concat(coloc.toList())
}