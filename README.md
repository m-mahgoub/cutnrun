### **!!!! This pipeline is still under development !!!!**
<br>

# cutnrun

cutnrun is a Nextflow pipeline for CUT&RUN analysis.

## Installation

```bash
git clone https://github.com/m-mahgoub/cutnrun.git
```

## Pipeline summary
1. Quality check (fastqc)
2. Generate bowtie2 index for reference genome (bowtie2 )
3. Mapping with bowtie2 (bowtie2 and samtools)
4. Calling peaks (macs2)
5. Generate bigwig files using (deeptools bamcoverage)
6. Generate heatmaps (deeptools computeMatrix and plotHeatmaps)
* The command-line options for these tools can be conveniently modified in nextflow.config file




## Usage
### Input files:
1. fastq raw reads (either paired-end or single-end), provided as a full path in the samplesheet.csv files
2. Path for reference genome fasta file provided as a parameter in nextflow.config file
3. Define the plotting strategy for the desired heatmaps in YAML format in “heatmap_blueprint.yaml” file (as shown in the sample file). The user can add as many plots as required, and specify:
     1) The name of the plot
     2) The paths and labels (optional) of the bed files for regions of plotting
     3) The paths and labels (optional) of the bigwig files to be plotted.
    * These files in ii and iii can be either a sample in the current pipeline (only ID is provided), files in the local path of the environment running the Nextflow pipeline, or remote files in servers.

### Requirements:
This pipeline requires the following dependencies to be installed and avaialable in the path of the environment running the pipeline:
1. **[Nextflow >= 20.07.1](https://www.nextflow.io/). Older Nextflow versions don't support DSL2 synatax used in this pipeline and they are not comaptible.**
2. [Singularity](https://sylabs.io/singularity) or [Docker](https://www.docker.com/).

### Running pipeline in local environment
with Singularity
```bash
nextflow run test.nf -profile singularity
```
or with Docker
```bash
nextflow run test.nf -profile docker
```

### Running pipeline in SLURM HPC cluster
with Singularity
```bash
nextflow run test.nf -profile slurm,singularity
```
or with Docker
```bash
nextflow run test.nf -profile slurm,docker
```

### Test run: 
1. Clone the repository
2. Run the previous commands in any environment satisfying the requirements.
3. Output from test run should be generated in a directory named "results". If the test run is successful, "results" directory should be identical to the directory named "testOut"

## Credit:
- Sample Datasets used in this pipeline are from nf-core chipseq repository (https://github.com/nf-core/chipseq)
- Original Sample Data Sheet: https://raw.githubusercontent.com/nf-core/cutnrun-test-datasets/chipseq/design.csv
