# Sifaka_assembly
This repository contains code to assemble and analyze sifaka exomic data.

For questions, contact:
Tim Webster
Arizona State University

## Goals
The main goals of this project are to explore:
1) the success of human exome capture probes in Strepsirrhine primates (with a Catarrhine comparison)
2) the success of a number of different assembly and variant calling techniques in reference and non-reference based settings

## Data
This project uses exome data generated (using human capture probes) from:
- 4 sifakas (mother-daughter pair, two unrelated males)
- 4 Indian rhesus macaques (2 males, 2 females, all unrelated)

## Setup
The assembly/analysis pipeline is written in Snakemake.  It handles running a variety of external tools, as well as a few internal programs (e.g., mapqs).  To run analyses, I strongly suggest using Anaconda (with bioconda added as a channel) to manage packages and programs.  When analyses are complete, I will provide an environment.yml file with exact versions of each program.  Until then, use the following commands as a suggestion. *I assume you are working with a Linux operating system.*

* First, install [Anaconda](https://www.continuum.io/downloads) or [Miniconda](https://conda.io/miniconda.html) for Python version 3.6 or higher.

* Next, add relevant channels to conda (e.g., bioconda) with the following commands:
	```
	conda config --add channels r
	conda config --add channels conda-forge
	conda config --add channels bioconda
	```
* Create the environments we'll be working in and install required packages with the commands:

	```
	conda create --name Sifaka_snakemake python=3.6 bbmap bedtools htslib bwa fastqc freebayes samblaster samtools ensembl-vep

	conda create --name sifaka_platypus platypus-variant
	```

 These two commands will create two separate environments, one (Sifaka_snakemake) that will be used for the main pipeline and most commands and a second for Platypus variant calling because Platypus requires a different version of Python.

 * Download [GATK](https://software.broadinstitute.org/gatk/download/) (we used version 3.7)

 * Build ```scripts/mapqs```

1. Check to see if Go is installed on your system: ```go version```
	* If not, download and install Go ([instructions here](https://golang.org/doc/install))
2. Get the required go packages: ```go get github.com/biogo/hts/...```
3. Change to the scripts directory and compile: ```go build mapqs.go```

You can check ```mapqs``` by changing to the `scripts` directory and entering the command:

	```
	./mapqs -h
	```
If you get something that looks like:

	```
	Usage of ./mapqs:
	  -help
		   display help
	  -infile string
		   input bam file
	  -outfile string
		   output file for stats
	  -threads int
		   number of threads to use (0 = auto)
	```

You're all set.

## Running the pipeline
* First, load the main anaconda environment:
```
source activate Sifaka_snakemake
```

* Change any relevant variables (namely paths) at the top of ```Snakefile```

* From the directory in which ```Snakefile``` resides, try the command:
```
snakemake -np
```
This should result in a long list of commands, wildcards, jobs to run, etc., in
green with a final summary of job counts per rule.  If not, check the directory you're
in and read any errors to sort out problems you might have.  For example: did you
update the variables in ```Snakefile``` correctly? Are your fastq files named
the same way that they are in the script (e.g., sampleid_read1.fastq.gz and so on)?

* If everything worked and you got a list of jobs to be run, you can run everything
with the command:
```
snakemake
```
Obviously, you might want to customize things a bit, by telling Snakemake to run
jobs in parallel, deploy on a cluster, etc.  You can find out much more information
in the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html).
