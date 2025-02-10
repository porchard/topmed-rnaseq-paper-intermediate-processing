#!/usr/bin/env nextflow

nextflow.enable.dsl=2
SUSIE_GLOB = params.susie_glob
CS_GLOB = params.cs_glob
CONTAINER = params.container
BIN = params.bin

process make_R_object {

    maxForks 50
    container "$CONTAINER"
    publishDir "${params.results}/coloc-in/${tissue}"
    memory '30 GB'
    time '2h'

    input:
    tuple val(tissue), path(susie), path(cs)

    output:
    path("${tissue}.*.rda")

    """
    cut -f1,5 $cs | uniq | sort | uniq | perl -pe 's/\t/\tL/' > cs.txt
    python3.8 ${BIN}/susie-pickle r-object --pickle $susie --keep cs.txt --prefix ${tissue}.
    """

}


workflow {
    susie_in = Channel.fromPath(SUSIE_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, susie
    cs_in = Channel.fromPath(CS_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // tissue, cs file
    make_R_object(susie_in.combine(cs_in, by: 0))
}