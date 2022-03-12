nextflow.enable.dsl=2

/* Channel.fromFilePairs(params.reads)
    .view() */


// println(params.reads)

process fastqc {
    publishDir params.outdir, mode:'copy'
    input:
        tuple val(sampleID), path(reads)
    output:
        // path "${sampleID}.txt"
        stdout emit: cmd

    script:
    """
    echo fastqc ${reads[0]} $sampleID
    echo fastqc ${reads[1]} $sampleID
    """
}

workflow {
    reads = Channel.fromFilePairs(params.reads)
    fastqc_ch = fastqc(reads)
    fastqc_ch.view()
}