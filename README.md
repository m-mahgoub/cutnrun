!!!! This pipeline is still under development !!!!


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
2. Path for reference genome fasta field provided as a parameter in nextflow.config file
3. Define the plotting strategy for the desired heatmaps in YAML format in “deeptools_user_metadata.yaml” file (as shown in the sample file). The user can add as many plots as required, and specify:
     1) The name of the plot
     2) The paths and labels (optional) of the bed files for regions of plotting
    3) The paths and labels (optional) of the bigwig files to be plotted.
    * These files in 2 and 3 can be either a sample in the current pipeline (only ID is provided), files in the local path of the environment running the Nextflow pipeline, or remote files in servers.

### Requirements:
It is required the following dependencies are installed in the environment running the pipeline:
1. Nextflow >= 20.07.1.
2. Conda and Singularity (or Docker).

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


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)