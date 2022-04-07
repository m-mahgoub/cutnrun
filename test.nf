nextflow.enable.dsl=2

if (params.single_end) { 
     Channel.fromPath(params.sample_sheet)
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true) ] ] }
    .set { raw_reads_fastq_ch }
} else { 
    Channel.fromPath(params.sample_sheet)
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ] }
    .set { raw_reads_fastq_ch }
                        
}


process get_controls {
    publishDir "$params.outdir/metadata", mode:'copy', pattern: '*'
    cpus 2
    memory '4 GB'
    input:
        path sample_sheet
    output:
        path 'chip_input_match.csv'
 
    script:
        """
        get_controls.py $sample_sheet chip_input_match.csv
        """

}


process get_deeptools_user_metadata {
    publishDir "$params.outdir/metadata", mode:'copy', pattern: '*'
    cpus 2
    memory '4 GB'
    input:
        path deeptools_user_metadata_yaml
    output:
        path 'bash_strings_for_deeptools.txt'
        path 'files_paths_for_deeptools.txt'
    script:
        """
        get_deeptools_user_metadata.py $deeptools_user_metadata_yaml bash_strings_for_deeptools.txt files_paths_for_deeptools.txt
        """

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
        tuple val(sampleID), val(controlID)
        path bamFile
    output:
        path "${sampleID}.bed"
 
    script:
    if( params.single_end )
        """
        macs2 callpeak -t ${sampleID}.sorted.bam -c ${controlID}.sorted.bam -n ${sampleID} -f "BAM" --nomodel --extsize 200 ${params.processes_options.macs2}
        mv ${sampleID}_peaks.narrowPeak ${sampleID}.bed
        """
    else
        """
        macs2 callpeak -t ${sampleID}.sorted.bam -c ${controlID}.sorted.bam -n ${sampleID} -f "BAMPE" ${params.processes_options.macs2}
        mv ${sampleID}_peaks.narrowPeak ${sampleID}.bed
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

process name {
    publishDir "$params.outdir/deeptools/matrix", mode:'copy', pattern: '*'
    cpus 2
    memory '4 GB'
    input:
        path X
        tuple val(X), path(X)
    output:
        path "X."
        path "X."
 
    script:
    if( params.X == XX )
        """
        echo X
        """
    else
        """
        echo X
        """
}



workflow {
    get_controls(Channel.fromPath(params.sample_sheet))
    
    // make channle that emmits pairs of sample, control
    get_controls.out[0]
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.sample_id, row.control] }
    .set { sample_control_pair_ch }

    get_deeptools_user_metadata(Channel.fromPath(params.deeptools_yaml))
    get_deeptools_user_metadata.out[0].view()

    // make channle that emmits parameters and external files paths for deeptools plotting
    get_deeptools_user_metadata.out[0]
    .splitCsv(header:true, sep:'\t')
    .map { row -> [ row.plot_name, [row.bed_files_inLine, row.bigwig_files_inLine, row.bed_labels_inLine_with_quotes, row.bigwig_labels_inLine_with_quotes ]] }
    .set { strings_for_deeptools_ch }
	// strings_for_deeptools_ch.view()

    // make channle that emmits files paths for deeptools
    get_deeptools_user_metadata.out[1]
    .splitCsv(header:true, sep:'\t')
    .map { row -> [ row.plot_name, row.paths_to_include] }
    .set { paths_for_deeptools_ch }
	// paths_for_deeptools_ch.view()
    paths_for_deeptools_ch.join(strings_for_deeptools_ch)
    .set { deeptools_meta_ch }
    deeptools_meta_ch.view()
    



    fastqc(raw_reads_fastq_ch)
    bowtie2_index(Channel.fromPath(params.genome))
    bowtie_mapping(bowtie2_index.out.first(), raw_reads_fastq_ch)
    callpeaks(sample_control_pair_ch,bowtie_mapping.out[0].collect())
    make_bigwig(bowtie_mapping.out[0], bowtie_mapping.out[1])
}