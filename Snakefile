import os

configfile: "sifaka_config.json"

fastq_directory = "/home/thwebste/Data/sifaka"
temp_directory = "temp/"

pcoq_1_path = "reference/GCF_000956105.1_Pcoq_1.0_genomic.fna"
pcoq_1_prefix = "reference/GCF_000956105.1_Pcoq_1.0_genomic"
hg38_path = "reference/Homo_sapiens_assembly38.fasta"
hg38_prefix = "reference/Homo_sapiens_assembly38"

# xyalign_path =
samtools_path = "samtools"
bwa_path = "bwa"
bbmerge_sh_path = "bbmerge.sh"
bbduksh_path = "bbduk.sh"
fastqc_path = "fastqc"
samblaster_path = "samblaster"
gatk_path = "/home/thwebste/Tools/GenomeAnalysisTK_37.jar"

all_samples = config["males"] + config["females"]

rule all:
	input:
		expand("fastqc/{fq_prefix}_fastqc.html", fq_prefix=config["fastq_prefixes"]),
		expand("adapters/{sample}.adapters.fa", sample=all_samples),
		expand("processed_bams/{sample}.hg38.sorted.mkdup.bam", sample=all_samples),
		expand("processed_bams/{sample}.pcoq.sorted.mkdup.bam", sample=all_samples),
		expand("stats/{sample}.pcoq.sorted.mkdup.bam.stats", sample=all_samples),
		expand("stats/{sample}.hg38.sorted.mkdup.bam.stats", sample=all_samples),
		expand("callable_sites/{sample}.hg38.ONLYcallablesites.bed", sample=all_samples),
		expand("callable_sites/{sample}.pcoq.ONLYcallablesites.bed", sample=all_samples)

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

rule prepare_reference_hg38:
	input:
		hg38_path
	output:
		fai = hg38_path + ".fai",
		amb = hg38_path + ".amb",
		dict = hg38_prefix + ".dict"
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
		fq1 = os.path.join(fastq_directory, "{sample}_read1.fastq.gz"),
		fq2 = os.path.join(fastq_directory, "{sample}_read2.fastq.gz")
	output:
		"adapters/{sample}.adapters.fa"
	params:
		bbmerge_sh = bbmerge_sh_path
	shell:
		"{params.bbmerge_sh} in1={input.fq1} in2={input.fq2} outa={output} reads=1m"

rule trim_adapters_paired_bbduk:
	input:
		fq1 = os.path.join(fastq_directory, "{sample}_read1.fastq.gz"),
		fq2 = os.path.join(fastq_directory, "{sample}_read2.fastq.gz")
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx3g in1={input.fq1} in2={input.fq2} out1={output.out_fq1} out2={output.out_fq2} ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tbe qtrim=rl trimq=10"

rule fastqc_analysis_trimmed_paired:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		html1 = "fastqc/{sample}_trimmed_read1_fastqc.html",
		html2 = "fastqc/{sample}_trimmed_read2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc {input.fq1} {input.fq2}"

rule map_and_process_trimmed_reads_hg38:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		fai = hg38_path + ".fai",
		ref = hg38_path
	output:
		"processed_bams/{sample}.hg38.sorted.mkdup.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samtools = samtools_path,
		samblaster = samblaster_path
	threads: 4
	shell:
		" {params.bwa} mem -t {threads} -R "
		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samblaster} "
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

rule map_and_process_trimmed_reads_pcoq:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		fai = pcoq_1_path + ".fai",
		ref = pcoq_1_path
	output:
		"processed_bams/{sample}.pcoq.sorted.mkdup.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samtools = samtools_path,
		samblaster = samblaster_path
	threads: 4
	shell:
		" {params.bwa} mem -t {threads} -R "
		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samblaster} "
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

rule bam_stats:
	input:
		"processed_bams/{sample}.{genome}.sorted.mkdup.bam"
	output:
		"stats/{sample}.{genome}.sorted.mkdup.bam.stats"
	shell:
		"samtools stats {input} | grep ^SN | cut -f 2- > {output}"

rule index_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.mkdup.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.mkdup.bam.bai"
	shell:
		"samtools index {input}"

rule generate_callable_sites_hg38:
	input:
		ref = hg38_path,
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai"
	output:
		"callable_sites/{sample}.{genome}.callablesites"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk,
		summary = "stats/{sample}.{genome}.callable.summary"
	shell:
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T CallableLoci -R {input.ref} -I {input.bam} --minDepth 4 --summary {params.summary} -o {output}"

rule extract_callable_sites:
	input:
		"callable_sites/{sample}.{genome}.callablesites"
	output:
		"callable_sites/{sample}.{genome}.ONLYcallablesites.bed"
	shell:
		"sed -e '/CALLABLE/!d' {input} > {output}"

rule generate_callable_sites_pcoq:
	input:
		ref = pcoq_1_path,
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai"
	output:
		"callable_sites/{sample}.{genome}.callablesites"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk,
		summary = "stats/{sample}.{genome}.callable.summary"
	shell:
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T CallableLoci -R {input.ref} -I {input.bam} --minDepth 4 --summary {params.summary} -o {output}"

rule extract_callable_sites:
	input:
		"callable_sites/{sample}.{genome}.callablesites"
	output:
		"callable_sites/{sample}.{genome}.ONLYcallablesites.bed"
	shell:
		"sed -e '/CALLABLE/!d' {input} > {output}"
