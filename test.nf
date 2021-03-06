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
    container  "mdivr/conda-nf-cutnrun:v0.1"
    time '10m'
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


process heatmap_blueprint {
    time '10m'
    container  "mdivr/conda-nf-cutnrun:v0.1"
    publishDir "$params.outdir/metadata", mode:'copy', pattern: '*'
    cpus 2
    memory '4 GB'
    input:
        path heatmap_blueprint_yaml
    output:
        path 'bash_strings_for_deeptools.txt'
        path 'files_paths_for_deeptools.txt'
    script:
        """
        get_heatmap_blueprint.py $heatmap_blueprint_yaml bash_strings_for_deeptools.txt files_paths_for_deeptools.txt
        """

}


process fastqc {
    time '30m'
    container  "quay.io/biocontainers/fastqc:0.11.9--0"
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
    time '30m'
    container  "quay.io/biocontainers/bowtie2:2.4.5--py38he5f0661_1"
    publishDir "$params.outdir/", mode:'copy'
    cpus 2
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
    time '30m'
    container  "alexeyebi/bowtie2_samtools:latest"
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
    time '30m'
    container  "quay.io/biocontainers/macs2:2.2.7.1--py38hbff2b2d_4"
    publishDir "$params.outdir/peaks/narrow", mode:'copy', pattern: '*.bed'
    cpus 2
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
    time '30m'
    container  "quay.io/biocontainers/deeptools:3.5.1--py_0"
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
        bamCoverage -p ${task.cpus} --bam $bamFile -o ${bamFile.getSimpleName()}.bigwig ${params.processes_options.bamcoverage}
        """

}

process plotHeatmap {
    time '30m'
    container  "quay.io/biocontainers/deeptools:3.5.1--py_0"
    publishDir "$params.outdir/deeptools/matrix", mode:'copy', pattern: '*.gz'
    publishDir "$params.outdir/deeptools/heatmaps", mode:'copy', pattern: '*.png'
    cpus 2
    memory '4 GB'
    input:
        path callpeaks_out
        path make_bigwig_out
        path (local_remote_files)
        tuple val(plot_name), val (bed_files_inLine), val (bigwig_files_inLine), val(bed_labels_inLine_with_quotes),  val(bigwig_labels_inLine_with_quotes)
    output:
        path "${plot_name}.gz"
        path "${plot_name}.png"
 
    script:
    """
    computeMatrix reference-point --referencePoint center -R $bed_files_inLine -S $bigwig_files_inLine -o ${plot_name}.gz -p ${task.cpus} ${params.processes_options.computeMatrix}
    plotHeatmap -m ${plot_name}.gz -out ${plot_name}.png ${params.processes_options.plotHeatmap} --regionsLabel $bed_labels_inLine_with_quotes --samplesLabel $bigwig_labels_inLine_with_quotes
    """

}



workflow {
    get_controls(Channel.fromPath(params.sample_sheet))
    
    // make channle that emmits pairs of sample, control
    get_controls.out[0]
    .splitCsv(header:true, sep:',')
    .map { row -> [ row.sample_id, row.control] }
    .set { sample_control_pair_ch }

    heatmap_blueprint(Channel.fromPath(params.heatmap_blueprint_yaml))


    // make channle that emmits deepttols beds and bigwig options as string
    heatmap_blueprint.out[0]
    .splitCsv(header:true, sep:'\t')
    .map { row -> [ row.plot_name, row.bed_files_inLine, row.bigwig_files_inLine, row.bed_labels_inLine_with_quotes, row.bigwig_labels_inLine_with_quotes] }
    .set { strings_for_deeptools_ch }
    // make channle that emmits local and remote files paths for deeptools
    heatmap_blueprint.out[1]
    .splitCsv(header:false, sep:'\t')
    .collect()
    .set { paths_for_deeptools_ch }

    fastqc(raw_reads_fastq_ch)
    bowtie2_index(Channel.fromPath(params.genome))
    bowtie_mapping(bowtie2_index.out.first(), raw_reads_fastq_ch)
    callpeaks(sample_control_pair_ch,bowtie_mapping.out[0].collect())
    make_bigwig(bowtie_mapping.out[0], bowtie_mapping.out[1])
    plotHeatmap(callpeaks.out.collect(), make_bigwig.out.collect(), paths_for_deeptools_ch, strings_for_deeptools_ch)
}
