import os

configfile: "sifaka_config.json"

fastq_directory = "~/Data/sifaka"

pcoq_1_path = "reference/GCF_000956105.1_Pcoq_1.0_genomic.fna"
pcoq_1_prefix = "reference/GCF_000956105.1_Pcoq_1.0_genomic

# xyalign_path = 
samtools_path = "samtools"
bwa_path = "bwa"
bbmerge_sh_path = "bbmerge.sh"
bbduksh_path = "bbduk.sh"
fastqc_path = "fastqc"

all_samples = config["males"] + config["females"]

rule all:
	input:
		expand("fastqc/{fq_prefix}_fastqc.html", fq_prefix=config["fastq_prefixes"]),
		expand("adapters/{sample}.adapters.fa", sample=all_samples)

rule prepare_reference_pcoq_1:
	input:
		pcoq_1_path
	output:
		fai = pcoq_1_path + ".fai",
		amb = pcoq_1_path + ".amb",
		dict = pcoq_1_prefix + ".dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell("{params.samtools} faidx {input}")
		# .dict
		shell("{params.samtools} dict -o {output.dict} {input}")
		# bwa
		shell("{params.bwa} index {input}")

rule fastqc_analysis:
	input:
		os.path.join(fastq_directory, "{fq_prefix}.fastq.gz")
	output:
		"fastqc/{fq_prefix}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc {input}"

rule adapter_discovery:
	input:
		fq1 = os.path.join(fastq_directory, "{sample}_read1.fastq.gz")
		fq2 = os.path.join(fastq_directory, "{sample}_read2.fastq.gz")
	output:
		"adapters/{sample}.adapters.fa"
	params:
		bbmerge_sh = bbmerge_sh_path
	shell:
		"{params.bbmerge_sh} in1={input.fq1} in2={input.fq2} outa={output} reads=1m"
