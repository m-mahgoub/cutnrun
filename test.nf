nextflow.enable.dsl=2

/* Channel.fromFilePairs(params.reads)
    .view() */


// println(params.reads)

// process fastqc {
//     publishDir "$params.outdir/fastqc", mode:'copy'
//     input:
//         tuple val(sampleID), path(reads)
//     output:
//         path "${sampleID}_1.tx"

//     script:
//     """
//     echo fastqc ${reads[0]} > ${sampleID}_1.tx
//     echo fastqc ${reads[1]} > ${sampleID}_2.tx
//     """
// }

process fastqc {
    publishDir "$params.outdir/fastqc", mode:'copy'
    cpus 1
    memory '4 GB'
    input:
        tuple val(sampleID), path(reads)
    output:
        path "${sampleID}_logs"

    script:
    """
    mkdir ${sampleID}_logs
    fastqc ${reads[0]} ${reads[1]} --outdir ${sampleID}_logs
    """
}

process bowtie2_index {
    publishDir "$params.outdir/bowtie_index", mode:'copy'
    cpus 2
    memory '4 GB'
    input:
        path(genome_fasta)
    output:
        path ("index.*")
    script:
    """
    bowtie2-build --threads ${task.cpus} $genome_fasta index
    """

}


// process bowtie_mapping {
//     publishDir "$params.outdir/bowtie_mapping", mode:'copy'
//     input:
//         tuple val(sampleID), path(reads)
//     output:
//         path "${sampleID}.bam"
//     script:
//     """
//     echo bowtie2
//     """
// }

workflow {
    reads = Channel.fromFilePairs(params.reads)
    fastqc_ch = fastqc(reads)
    genome_ch = Channel.fromPath(params.genome)
    bowtie2_index_ch = bowtie2_index(genome_ch)
    // bowtie_mapping_ch = bowtie_mapping(reads)
}

