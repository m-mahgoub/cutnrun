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
    memory '4 GB'
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
    memory '4 GB'
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
    publishDir "$params.outdir/bowtie_mapping/bams", mode:'copy', pattern: '*.ba*'
    publishDir "$params.outdir/bowtie_mapping/logs", mode:'copy', pattern: '*.log'
    cpus 2
    memory '4 GB'
    input:
        path index
        tuple val(sampleID), path(reads)
    output:
        path "${sampleID}.sorted.bam"
        path "${sampleID}.sorted.bai"
        path "${sampleID}.log"

    script:
    if( params.single_end )
        """
        bowtie2 -x $index/index -U ${reads[0]}  --threads ${task.cpus} ${params.processes_options.bowtie2} > ${sampleID}.sam 2>${sampleID}.log
        samtools view -u ${sampleID}.sam | samtools sort > ${sampleID}.sorted.bam
        samtools index ${sampleID}.sorted.bam ${sampleID}.sorted.bai
        """
    else
        """
        bowtie2 -x $index/index -1 ${reads[0]}  -2 ${reads[1]} --threads ${task.cpus} ${params.processes_options.bowtie2} > ${sampleID}.sam 2>${sampleID}.log
        samtools view -u ${sampleID}.sam | samtools sort > ${sampleID}.sorted.bam
        samtools index ${sampleID}.sorted.bam ${sampleID}.sorted.bai
        """

}

process callpeaks {
    publishDir "$params.outdir/peaks/narrow", mode:'copy', pattern: '*.bed'
    cpus 4
    memory '4 GB'
    input:
        path bamFile
    output:
        path "${bamFile.getSimpleName()}.bed"
 
    script:
    if( params.single_end )
        """
        macs2 callpeak -t $bamFile -n ${bamFile.getSimpleName()} -f "BAM" --nomodel --extsize 200 ${params.processes_options.macs2}
        mv ${bamFile.getSimpleName()}_peaks.narrowPeak ${bamFile.getSimpleName()}.bed
        """
    else
        """
        macs2 callpeak -t $bamFile -n ${bamFile.getSimpleName()} -f "BAMPE" ${params.processes_options.macs2}
        mv ${bamFile.getSimpleName()}_peaks.narrowPeak ${bamFile.getSimpleName()}.bed
        """
}

process make_bigwig {
    publishDir "$params.outdir/deeptools/bigwig", mode:'copy', pattern: '*.bigwig'
    cpus 2
    memory '4 GB'
    input:
        path bamFile
        path bamIndexFile
    output:
        path "${bamFile.getSimpleName()}.bigwig"

    script:
        """
        bamCoverage --bam $bamFile -o ${bamFile.getSimpleName()}.bigwig ${params.processes_options.bamcoverage}
        """

}


workflow {
    fastqc(raw_reads_fastq_ch)
    bowtie2_index(Channel.fromPath(params.genome))
    bowtie_mapping(bowtie2_index.out, raw_reads_fastq_ch)
    callpeaks(bowtie_mapping.out[0])
    make_bigwig(bowtie_mapping.out[0], bowtie_mapping.out[1])
}