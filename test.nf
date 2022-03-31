nextflow.enable.dsl=2

if (params.single_end) { 
    raw_reads_fastq_ch = Channel.fromPath(params.sample_sheet)
                        .splitCsv(header:true, sep:',')
                        .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true) ] ] }
} else { 
    raw_reads_fastq_ch = Channel.fromPath(params.sample_sheet)
                        .splitCsv(header:true, sep:',')
                        .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ] }
                        
}

process fastqc {
    publishDir "$params.outdir/fastqc", mode:'copy'
    cpus 2
    memory '8 GB'
    input:
        tuple val(sampleID), path(reads)
    output:
        path "${sampleID}_logs"


    script:
    if( params.single_end )
        """
        mkdir ${sampleID}_logs
        fastqc ${reads[0]} --outdir ${sampleID}_logs
        """
    else
        """
        mkdir ${sampleID}_logs
        fastqc ${reads[0]} ${reads[1]} --outdir ${sampleID}_logs
        """
}

process bowtie2_index {
    publishDir "$params.outdir/bowtie_index", mode:'copy'
    cpus 4
    memory '8 GB'
    input:
        path(genome_fasta)
    output:
        path ("bowtie_index")
    script:
    """
    bowtie2-build --threads ${task.cpus} $genome_fasta index
    mkdir bowtie_index
    mv index.* bowtie_index
    """

}


process bowtie_mapping {
    publishDir "$params.outdir/bowtie_mapping/bams", mode:'copy', pattern: '*.bam'
    publishDir "$params.outdir/bowtie_mapping/logs", mode:'copy', pattern: '*.log'
    cpus 2
    memory '4 GB'
    input:
        path index
        tuple val(sampleID), path(reads)
    output:
        path "${sampleID}.sorted.bam"
        path "${sampleID}.log"
    script:
    """
    bowtie2 -x $index/index -1 ${reads[0]}  -2 ${reads[1]} --threads ${task.cpus} ${params.processes_options.bowtie2} > ${sampleID}.sam 2>${sampleID}.log
    samtools view -u ${sampleID}.sam | samtools sort > ${sampleID}.sorted.bam

    """
    script:
    if( params.single_end )
        """
        bowtie2 -x $index/index -U ${reads[0]}  --threads ${task.cpus} ${params.processes_options.bowtie2} > ${sampleID}.sam 2>${sampleID}.log
        samtools view -u ${sampleID}.sam | samtools sort > ${sampleID}.sorted.bam
        """
    else
        """
        bowtie2 -x $index/index -1 ${reads[0]}  -2 ${reads[1]} --threads ${task.cpus} ${params.processes_options.bowtie2} > ${sampleID}.sam 2>${sampleID}.log
        samtools view -u ${sampleID}.sam | samtools sort > ${sampleID}.sorted.bam
        """

}

workflow {
    fastqc_ch = fastqc(raw_reads_fastq_ch)
    genome_ch = Channel.fromPath(params.genome)
    bowtie2_index(genome_ch)
    bowtie_mapping_ch = bowtie_mapping(bowtie2_index.out, raw_reads_fastq_ch)
    // ch_raw_reads_fastq.view()
    
}


