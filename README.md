# Sifaka_assembly
This repository contains code to assemble and analyze sifaka exome data.

For questions, contact:
Tim Webster
Department of Anthropology, University of Utah

## Citation
Please cite our preprint:

Webster TH, Guevara EE, Lawler RR, Bradley BJ. 2018. Successful exome capture and sequencing in lemurs using human baits.
bioRxiv 490839. doi: https://doi.org/10.1101/490839

## Goals
The main goal of this project is to explore the success of human exome capture probes in strepsirrhine primates (with a catarrhine comparison)

## Data
This project uses exome data generated (using human capture probes) from:

*Primary dataset*
- 4 sifakas (mother-daughter pair, two unrelated males; dataset Sifaka1 in the paper): M288, M418, F249, F406
- 4 Indian rhesus macaques (2 males, 2 females, all unrelated; dataset Macaque1 in the paper): WI055, WI059, WI056, WI057

*Secondary dataset*
- 4 additional sifakas (2 males, 2 females, all unrelated; dataset Sifaka2 in the paper): 137, 161, 184, 252

We used the primary dataset to directly compare exome capture success between rhesus macaques and sifakas. Human exome probes are known to work well in macaques, so they provide an ideal comparison.

We used the four additional sifakas in the secondary dataset to understand how exome bait design might affect success. The primary dataset was captured using NimbleGen SeqCap EZ Exome v2 probes, while the secondary dataset was captured using IDT xGen Exome Research Panel probes.

*Reference genomes*
We downloaded reference genomes to the ``reference`` directory.

We obtained pcoq1 from https://www.ncbi.nlm.nih.gov/genome/24390

We obtained mmul from https://www.ncbi.nlm.nih.gov/genome/215?genome_assembly_id=259055

We obtained hg38 from the Broad Institute's resource bundle: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0?pli=1

## Setup
The assembly/analysis pipeline is written in Snakemake.  It handles running a variety of external tools, as well as a few internal programs (e.g., mapqs).  To run analyses, I strongly suggest using Anaconda (with bioconda added as a channel) to manage packages and programs.  When analyses are complete, I will provide an environment.yml file with exact versions of each program.  Until then, use the following commands as a suggestion. *I assume you are working with a Linux operating system.*

* First, install [Anaconda](https://www.continuum.io/downloads) or [Miniconda](https://conda.io/miniconda.html) for Python version 3.6 or higher.

* Next, add relevant channels to conda (e.g., bioconda) with the following commands:
	```
	conda config --add channels r
	conda config --add channels defaults
	conda config --add channels conda-forge
	conda config --add channels bioconda
	```
* Create the environments we'll be working in and install required packages with the commands:

	```
	conda create --name Sifaka_snakemake python=3.6 bbmap bedtools htslib bwa fastqc freebayes samblaster samtools bedops cyvcf2 ensembl-vep

	```

 This command will create an environment (Sifaka_snakemake) that will contain all programs (except mapqs - see below) required for the pipeline.

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
* The pipeline is designed to work with the original fastq file names, not those used by SRA. So, if you're planning on using the fastq files downloaded from SRA make sure filenames are changed to:

| SRA ID | Expected filenames |
| ------ | ------------------ |
| SRR6277165 | F249_read1.fastq.gz, F249_read2.fastq.gz |
| SRR6277166 | F406_read1.fastq.gz, F406_read2.fastq.gz |
| SRR6277167 | M288_read1.fastq.gz, M288_read2.fastq.gz |
| SRR6277168 | M418_read1.fastq.gz, M418_read2.fastq.gz |
| SRR7219233 | 137_read1.fastq.gz, 137_read2.fastq.gz |
| SRR7219232 | 161_read1.fastq.gz, 161_read2.fastq.gz |
| SRR7219231 | 184_read1.fastq.gz, 184_read2.fastq.gz |
| SRR7219230 | 252_read1.fastq.gz, 252_read2.fastq.gz |
| SRR5085658 | WI055_read1.fastq.gz, WI055_read2.fastq.gz |
| SRR5085657 | WI056_read1.fastq.gz, WI056_read2.fastq.gz |
| SRR5085688 | WI057_read1.fastq.gz, WI057_read2.fastq.gz |
| SRR5085648 | WI059_read1.fastq.gz, WI059_read2.fastq.gz |

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

## Steps in the pipline
The following text walks through and describes all of the steps in the pipeline. Text like ``rule <something>`` is directly referring to a specific rule in the Snakefile.

### Step 1: Preliminary setup (before Snakemake pipeline)
Before running the pipeline, we prepared conda environments (see above) and downloaded reference genomes. For the reference genomes, we used Mmul8.0.1 (rhesus macaque; downloaded from NCBI), Pcoq1.0 (Coquerel's sifaka; dowloaded from NCBI), and hg38 (human; downloaded from GATK's resource bundle).

Note also that the config file contains data for rhemac2, a previous version of the rhesus reference genome that we used in early stages of analysis. However, as we did not use this genome in our final analyses, we removed it from ``rule all`` in ``Snakefile``.

### Step 2: Prepare reference genomes
In the first rule, ``rule prepare_reference``, we first create a symbolic link, which provides a naming pattern that's much easier to work with in Snakemake. We write everything to a new directory ``new_reference`` to ensure that there aren't name clashes with the original reference genomes.

This rule then creates fasta indices (.fai), sequence dictionaries (.dict), and BWA indices (variety of extensions) that are needed for tools downstream.

### Step 3: Download, process, and prepare genome annotations
In this series of rules, we first download reference annotations from NCBI in ``rule get_annotation``. Note that this will not work for rhemac2 as written (which won't matter if you're only interested in the final version of the pipeline, which uses Mmul 8.0.1 instead).

We then either extract features from the gff file, or infer them via ``bedtools subtract``.

Finally, we convert hg38 scaffold names to those used in the reference we downloaded from GATK in ``rule convert_gff_coordinates``.

The final output of these rules is a series of gff files, each containing one type of genomic feature (exons, genes, cds, UTRs). Note that we will use ``bedtools subtract`` later in the pipeline to subtract genes, leaving intergenic regions.

### Step 4: Preliminary fastq quality assessment
In this step, in ``rule fastqc_analysis`` and ``rule multiqc_analysis``, we use fastqc and multiqc to calculate and summarize various quality metrics from all of the fastq files. The resulting multiqc report summarizing all of this information will be ``multiqc/multiqc_report.html``.

### Step 5: Adapter discovery
In early versions of the pipeline, we used ``rule adapter_discovery`` to identify adapters in the raw sequencing reads using BBMerge, which is part of the BBTools suite. All adapters were Illumina Universal adapters, which we include in ``misc/adapter_sequence.fa``. Note that running the full pipeline will use this rule and create the output, but it won't do anything with the output.

### Step 6: Fastq trimming and quality control
In ``rule trim_adapters_paired_bbduk``, we use BBDuk to trim the universal adapters (see Step 5 above) as well as apply some gentle quality trimming. The trimming parameters are based on, but slightly modified from BBDuk standard recommendations. The output of this rule is a pair of trimmed, gzipped fastq files located in the ``trimmed_fastqs`` directory.

### Step 7: Quality assessment of trimmed reads
We then repeat the analyses undertaken in Step 4 (above) on the trimmed fastqs (output of Step 6). This occurs, however, in new rules: ``rule fastqc_analysis_trimmed_paired`` and ``rule trimmed_multiqc_analysis``. The output is ``multiqc_trimmed/multiqc_report.html``.

### Step 8: Mapping reads, removing duplicates, and sorting BAM files
In this step, we take the trimmed reads output from Step 6, map them to their respective reference genome (human for macaques and sifakas, macaque for the macaques, and sifaka for the sifakas--see Steps 1 and 2), remove duplicates, fix read pair information, and sort the alignments. We do this all in ``rule map_and_process_trimmed_reads``.

For mapping, we use ``bwa mem`` and include read group information that we provided in the config file (read group ID, sample name, library, sequencing flowcell/lane, and sequencing platform). We then pipe the output to SAMBLASTER to mark, but not remove, duplicate reads. We then pipe that output to ``samtools fixmate`` to correct read pairing information and output in BAM format. Finally, we pipe to ``samtools sort`` to sort the reads in coordinate order and output a BAM file in the ``processed_bams`` directory.

### Step 9: Downsample BAM files
Our sequencing datasets differed quite a bit in the number of raw read generated. So, we downsampled our BAM files to contain the same number of nonduplicate reads.

To do this, we initially counted the number of nonduplicate reads (not in this pipeline) in each file. We then calculated the fraction of this number that would result in 50 million reads (important for the samtools command--see next paragraph). We include these values in the config file.

In ``rule downsample_bams``, we used ``samtools view`` to subset approximately 50 million nonuplicate reads from each sample. The ``-F 1024`` flag excludes duplicates, and the ``-s 0.(params.downsample_fraction)`` provides the fraction to include in the output. For example, if params.downsample_fraction was equal to 8942, then samtools would downsample to 89.42% of the total number of nonduplicate reads for that sample.

The output is in ``processed_bams``. For the rest of the pipeline, if ``downsampled`` is in a file name, it is using the downsampled data, while if ``unsampled`` is in a file name, it is containing data that has not been downsampled.

We further index all BAM files, both unsampled and downsampled, in ``rule index_bam`` with ``samtools index``.

### Step 10: BAM statistics
In ``rule bam_stats`` and ``rule bam_stats_on_passing_reads``, we use samtools to calculate a variety of metrics (raw reads, passing reads, duplicates, pairing, etc.) on all BAM files. All output is deposited in the ``stats`` directory.

### Step 11: Mean MAPQ
In this project, we're interested in the success of exome capture with human baits across the order primates. Previous studies have suggested that there's a major dropoff in success outside of Catarrhine primates. However, we hypothesized that mapping might actually be driving much of this pattern. Therefore, we needed to calculate mapping quality across BAM files.

We do this in ``rule mapq_check``, which runs a custom program written by Tim Webster to calculate mean and SD of MAPQ quickly across a BAM file. This program is written in GO and included in the distribution of this pipeline. See the Setup section above for information about compiling this program.

The output of this program is written to the ``stats`` directory (where Step 10 output is written as well). All output from this rule ends in ``.mapq.stats``. This output is input for plotting and analysis (see below).

### Step 12: Create bedgraph containing coverage across the genome
In this step, we use a series of rules to calculate depth of coverage (nonduplicate reads; minimum MAPQ of 20) at every site in the genome.

First, we need to create a "genome file" for bedtools. In ``rule make_bedtools_genome_file``, we use awk to select the first two columns of our .fai genome index files created in Step 2. These columns contain the scaffold name and scaffold length, respectively.

Next, we use this genome file and the BAM files created in Step 9 and Step 10 as input for ``rule genome_cov``, where we use the ``bedtools genomecov`` command to output a bedgraph of coverages per site. In this command, we first use ``samtools view`` to remove duplicate reads (``-F 1024``) and only include reads with a MAPQ of 20 or greater (``-q 20``).

Finally, we sort each bedgraph file using sort-bed, a utility in the BEDOPS suite.

All output is deposited in ``results``.

### Step 13: Compute histograms of coverage per region
The goal of this step is to output histograms of coverage in various regions of the genome.

In ``rule gff2bed``, we convert the genome annotation GFF files output in Step 3 to BED files using gff2bed, part of the BEDOPS suite. We use ``bedtools merge`` to ensure that regions aren't included multiple times (which commonly is the case in NCBI GFF files).

We then use bedtools to extract the regions from the bedgraphs output in Step 12. For most regions, this involves simply intersecting the bedgraph file with the region's file in ``rule bedtools_intersect_regions``. However, to get intergenic regions, we use bedtools to subtract genes from the bedgraph in ``rule bedtools_find_intergenic``. Note that in both cases we use ``bedtools merge`` as well to ensure that regions aren't included multiple times (which commonly is the case in NCBI GFF files).

Finally, we use a custom script, ``Compute_histogram_from_bed.py`` to calculate histograms from the extracted regions of the bedgraph files. The output of this step, in ``rule compute_histogram_from_bed`` and ``rule compute_histograph_from_bed_intergenic``, ends in ``.hist`` and is written to the ``results`` directory. These histograms are input for plotting and analysis (see below).

### Step 14: Calculate coverage breadth statistics for CDS
In ``rule cds_target_coverage``, instead of exploring depth of coverage statistics (like above), we use ``bedtools coverage`` to calculate, among other things, the fraction of each coding sequence (CDS) with at least one read mapped to it. For input, we use the output of ``rule gff2bed`` for CDS (Step 13) and the sorted bedgraph output from Step 12.

### Step 15: Find minimally callable sites
To speed up downstream variant calling, we first used GATK's CallableLoci tool to output a BED file containing, for each sample, sites classified by callability. This happens in ``rule generate_callable_sites``. We then use sed to extract callable regions from that BED file in ``rule extract_callable_sites``.

Then, in a series of rules beginning with "combine_callable_sites_" in their titles, we combine and merge regions into one BED file per species per reference genome containing every site that was minimally callable in at least one sample.

### Step 16: Generate GVCF files with GATK
Next, we use GATK's HaplotypeCaller to call variants for each sample, for each reference genome. We use the BED files containing minimally callable sites (output of Step 14) to constrain variant calling (for speed), and output a GVCF file that contains information about both variant and invariant sites.

### Step 17: Genotype GVCF files with GATK
Next, we use GATK's GenotypeGVCFs tool to genotype GVCF files for all samples from a given species, per reference genome. This outputs a raw VCF file per species, per reference genome. We use the flag ``--includeNonVariantSites`` to output records for every callable site in the genome.

### Step 18: Freebayes Variant Calling
In addition to GATK, we also used Freebayes to call variants. Again, we use the BED files output in Step 14 to restrict variant calling to regions considered minimally callable (for speed/efficiency). For Freebayes, we used default parameters and only output variant sites.

### Step 19: Zip and index VCF files
In ``rule zip_vcf`` we use bgzip to compress all VCF output from Steps 16 and 17. Next, in ``rule index_zipped_vcf``, we use tabix to index the zipped vcfs files.

### Step 20: Filter VCFs
We use a custom script, ``Filter_vcf.py``, to filter the compressed and indexed VCF files in ``rule filter_vcfs``. This script is built around the Python package cyvcf2, and we filtered for a minimum site quality of 30, sample depth of 8, and genotype quality of 30. We further required an allele to have at least three reads supporting it, and we only included sites callable in all samples of a given species/genome.

### Step 21: Zip and index filtered VCF files
Like Step 18, we use bgzip to compress and tabix to index the filtered vcf files.

### Step 22: Download HG38 Cache for VEP
In ``rule download_hg38_cache``, we download Ensembl's precompiled annotation database for hg38. We'll use gff files for the other two genomes.

### Step 23: Annotate VCFs
In ``rule vep_annotation`` we use Ensembl's Variant Effect Predictor (VEP) to annotate the VCF files output from Step 20. For hg38, we use the cache downloaded in Step 21, while we use GFFs and FASTA files for the other two genomes (Steps 2 and 3).

## Manual steps: Plotting
After the pipeline completed, we used R to manually plot and analyze the results in RStudio.

### Plot MAPQ results
We used the script ``mapq_figure_Rscript.R`` in the ``scripts`` directory to plot MAPQ results. It requires the packages ``cowplot`` and ``dpylr``. It takes as input a text file that we manually constructed as a part of our quality control checks. You can see the file we used in the ``misc`` directory. It's called ``mapping_results_all_samples.txt``.

### Plot coverage
We used ``plot_coverage_Rscript.R`` in the ``scripts`` directory for all coverage plots and analyses. It requires ``cowplot`` and ``dplyr``. This script needs to have the directory containing all of the ``.hist`` output files set as its working directory (if following the pipeline, it should be ``results``).

### Compare CDS coverage results
The ``cds_coverage_Rscript.R`` script in the ``scripts`` directory requires as input the files ending in ``.cdscoverage.bed`` and the working directory should be set to the directory containing these files (if following the pipeline, it should be ``results``). There are no required libraries if only calculating stats.
