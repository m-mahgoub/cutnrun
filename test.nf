nextflow.enable.dsl=2

/* Channel.fromFilePairs(params.reads)
    .view() */


// println(params.reads)

process fastqc {
    publishDir "$params.outdir/fastqc", mode:'copy'
    input:
        tuple val(sampleID), path(reads)
    output:
        path "${sampleID}_1.tx"

    script:
    """
    echo fastqc ${reads[0]} > ${sampleID}_1.tx
    echo fastqc ${reads[1]} > ${sampleID}_2.tx
    """
}

workflow {
    reads = Channel.fromFilePairs(params.reads)
    fastqc_ch = fastqc(reads)
}