/*
 * Define the default parameters 
 */

// nextflow.enable.dsl=2

// params.genome     = "$baseDir/data/genome.fa" 
// params.variants   = "$baseDir/data/known_variants.vcf.gz"
// params.blacklist  = "$baseDir/data/blacklist.bed"
// params.reads      = "$baseDir/data/reads/ENCSR000COQ1_{1,2}.fastq.gz" 
// params.results    = "results" 
// params.gatk       = "/opt/broad/GenomeAnalysisTK.jar" 


// /*
//  *  Parse the input parameters
//  */

// reads_ch        = Channel.fromFilePairs(params.reads)
// GATK            = params.gatk