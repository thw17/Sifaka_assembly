import os

configfile: "sifaka_config.json"

fastq_directory = "/home/thwebste/Data/sifaka"
temp_directory = "temp/"

pcoq_1_path = "reference/GCF_000956105.1_Pcoq_1.0_genomic.fna"
pcoq_1_prefix = "reference/GCF_000956105.1_Pcoq_1.0_genomic"
hg38_path = "reference/Homo_sapiens_assembly38.fasta"
hg38_prefix = "reference/Homo_sapiens_assembly38"

platypus_environment = "sifaka_platypus" # conda environment - platypus needs python 2.7

bbduksh_path = "bbduk.sh"
bbmerge_sh_path = "bbmerge.sh"
bedtools_path = "bedtools"
bgzip_path = "bgzip"
bwa_path = "bwa"
fastqc_path = "fastqc"
freebayes_path = "freebayes"
gatk = "/home/thwebste/Tools/GenomeAnalysisTK_37.jar"
platypus_path = "platypus"
samblaster_path = "samblaster"
samtools_path = "samtools"
tabix_path = "tabix"
# xyalign_path =

all_samples = config["males"] + config["females"]

rule all:
	input:
		expand(
			"fastqc/{fq_prefix}_fastqc.html", fq_prefix=config["fastq_prefixes"]),
		expand(
			"adapters/{sample}.adapters.fa", sample=all_samples),
		expand(
			"processed_bams/{sample}.hg38.sorted.mkdup.bam", sample=all_samples),
		expand(
			"processed_bams/{sample}.pcoq.sorted.mkdup.bam", sample=all_samples),
		expand(
			"stats/{sample}.pcoq.sorted.mkdup.bam.stats", sample=all_samples),
		expand(
			"stats/{sample}.hg38.sorted.mkdup.bam.stats", sample=all_samples),
		expand(
			"callable_sites/{sample}.hg38.ONLYcallablesites.bed", sample=all_samples),
		expand(
			"callable_sites/{sample}.pcoq.ONLYcallablesites.bed", sample=all_samples),
		expand(
			"stats/{sample}.pcoq.mapq.stats", sample=all_samples),
		expand(
			"stats/{sample}.hg38.mapq.stats", sample=all_samples),
		expand(
			"vcf/sifakas.pcoq.{caller}.raw.vcf.gz.tbi",
			caller=["freebayes", "platypus", "gatk"]),
		expand(
			"vcf/sifakas.hg38.{caller}.raw.vcf.gz.tbi",
			caller=["freebayes", "platypus", "gatk"])


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
		bam = "processed_bams/{sample}.hg38.sorted.mkdup.bam",
		bai = "processed_bams/{sample}.hg38.sorted.mkdup.bam.bai"
	output:
		"callable_sites/{sample}.hg38.callablesites"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk,
		summary = "stats/{sample}.hg38.callable.summary"
	shell:
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T CallableLoci -R {input.ref} -I {input.bam} --minDepth 4 --summary {params.summary} -o {output}"

rule generate_callable_sites_pcoq:
	input:
		ref = pcoq_1_path,
		bam = "processed_bams/{sample}.pcoq.sorted.mkdup.bam",
		bai = "processed_bams/{sample}.pcoq.sorted.mkdup.bam.bai"
	output:
		"callable_sites/{sample}.pcoq.callablesites"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk,
		summary = "stats/{sample}.pcoq.callable.summary"
	shell:
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T CallableLoci -R {input.ref} -I {input.bam} --minDepth 4 --summary {params.summary} -o {output}"

rule extract_callable_sites:
	input:
		"callable_sites/{sample}.{genome}.callablesites"
	output:
		"callable_sites/{sample}.{genome}.ONLYcallablesites.bed"
	shell:
		"sed -e '/CALLABLE/!d' {input} > {output}"

rule combine_callable_sites_pcoq:
	input:
		expand("callable_sites/{sample}.pcoq.ONLYcallablesites.bed", sample=all_samples)
	output:
		"callable_sites/combined.pcoq.ONLYcallablesites.bed"
	params:
		bedtools = bedtools_path
	shell:
		"cat {input} | sort -k1,1 -k2,2n | {params.bedtools} merge -i stdin > {output}"

rule combine_callable_sites_hg38:
	input:
		expand("callable_sites/{sample}.hg38.ONLYcallablesites.bed", sample=all_samples)
	output:
		"callable_sites/combined.hg38.ONLYcallablesites.bed"
	params:
		bedtools = bedtools_path
	shell:
		"cat {input} | sort -k1,1 -k2,2n | {params.bedtools} merge -i stdin > {output}"

rule mapq_check:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.mkdup.bam.bai"
	output:
		"stats/{sample}.{genome}.mapq.stats"
	threads: 4
	shell:
		"scripts/mapqs -infile {input.bam} -outfile {output} -threads {threads}"

rule gatk_gvcf_pcoq:
	input:
		ref = pcoq_1_path,
		bam = "processed_bams/{sample}.pcoq.sorted.mkdup.bam",
		bai = "processed_bams/{sample}.pcoq.sorted.mkdup.bam.bai",
		callable = "callable_sites/combined.pcoq.ONLYcallablesites.bed"
	output:
		"vcf/{sample}.pcoq.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads: 4
	shell:
		"java -Xmx16g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.callable} --emitRefConfidence GVCF -o {output}"

rule gatk_gvcf_hg38:
	input:
		ref = hg38_path,
		bam = "processed_bams/{sample}.hg38.sorted.mkdup.bam",
		bai = "processed_bams/{sample}.hg38.sorted.mkdup.bam.bai",
		callable = "callable_sites/combined.hg38.ONLYcallablesites.bed"
	output:
		"vcf/{sample}.hg38.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads: 4
	shell:
		"java -Xmx16g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.callable} --emitRefConfidence GVCF -o {output}"

rule genotype_gvcfs_pcoq:
	input:
		ref = pcoq_1_path,
		gvcfs = expand("vcf/{sample}.pcoq.g.vcf.gz", sample=all_samples)
	output:
		v = "vcf/sifakas.pcoq.gatk.raw.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads: 4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell("java -Xmx16g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T GenotypeGVCFs -R {input.ref} {variant_files} -o {output.v} --includeNonVariantSites")

rule genotype_gvcfs_hg38:
	input:
		ref = hg38_path,
		gvcfs = expand("vcf/{sample}.hg38.g.vcf.gz", sample=all_samples)
	output:
		v = "vcf/sifakas.hg38.gatk.raw.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads: 4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell("java -Xmx16g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T GenotypeGVCFs -R {input.ref} {variant_files} -o {output.v} --includeNonVariantSites")

rule platypus_variant_calling_pcoq:
	input:
		ref = pcoq_1_path,
		bams = expand("processed_bams/{sample}.pcoq.sorted.mkdup.bam", sample=all_samples),
		bais = expand("processed_bams/{sample}.pcoq.sorted.mkdup.bam.bai", sample=all_samples),
		callable = "callable_sites/combined.pcoq.ONLYcallablesites.bed"
	output:
		vcf = "vcf/sifakas.pcoq.platypus.raw.vcf"
	params:
		platypus = platypus_path,
		logfile = "logfiles/sifakas.pcoq.logfile"
	threads: 4
	run:
		bam_files = ",".join(input.bams)
		shell("source activate sifaka_platypus && {params.platypus} callVariants --output {output.vcf} --logFileName {params.logfile} --refFile {input.ref} --bamFiles {bam_files} --nCPU {threads} --targets {input.callable} --assemble=1 --outputRefCalls=1")

rule platypus_variant_calling_hg38:
	input:
		ref = hg38_path,
		bams = expand("processed_bams/{sample}.hg38.sorted.mkdup.bam", sample=all_samples),
		bais = expand("processed_bams/{sample}.hg38.sorted.mkdup.bam.bai", sample=all_samples),
		callable = "callable_sites/combined.hg38.ONLYcallablesites.bed"
	output:
		vcf = "vcf/sifakas.hg38.platypus.raw.vcf"
	params:
		platypus = platypus_path,
		logfile = "logfiles/sifakas.hg38.logfile"
	threads: 4
	run:
		bam_files = ",".join(input.bams)
		shell("source activate sifaka_platypus && {params.platypus} callVariants --output {output.vcf} --logFileName {params.logfile} --refFile {input.ref} --bamFiles {bam_files} --nCPU {threads} --targets {input.callable} --assemble=1 --outputRefCalls=1")

rule freebayes_pcoq:
	input:
		ref = pcoq_1_path,
		bams = expand("processed_bams/{sample}.pcoq.sorted.mkdup.bam", sample=all_samples),
		bais = expand("processed_bams/{sample}.pcoq.sorted.mkdup.bam.bai", sample=all_samples),
		callable = "callable_sites/combined.pcoq.ONLYcallablesites.bed"
	output:
		vcf = "vcf/sifakas.pcoq.freebayes.raw.vcf"
	params:
		freebayes = freebayes_path
	shell:
		"{params.freebayes} -f {input.ref} -v {output} --targets {input.callable} --report-monomorphic {input.bams}"

rule freebayes_hg38:
	input:
		ref = hg38_path,
		bams = expand("processed_bams/{sample}.hg38.sorted.mkdup.bam", sample=all_samples),
		bais = expand("processed_bams/{sample}.hg38.sorted.mkdup.bam.bai", sample=all_samples),
		callable = "callable_sites/combined.hg38.ONLYcallablesites.bed"
	output:
		vcf = "vcf/sifakas.hg38.freebayes.raw.vcf"
	params:
		freebayes = freebayes_path
	shell:
		"{params.freebayes} -f {input.ref} -v {output} --targets {input.callable} --report-monomorphic {input.bams}"

rule zip_vcf:
	input:
		vcf = "vcf/sifakas.{genome}.{caller}.raw.vcf"
	output:
		"vcf/sifakas.{genome}.{caller}.raw.vcf.gz"
	params:
		bgzip = bgzip_path
	shell:
		"{params.bgzip} {input.vcf}"

rule index_zipped_vcf:
	input:
		vcf = "vcf/sifakas.{genome}.{caller}.raw.vcf.gz"
	output:
		"vcf/sifakas.{genome}.{caller}.raw.vcf.gz.tbi"
	params:
		tabix = tabix_path
	shell:
		"{params.tabix} -p vcf {input}"
