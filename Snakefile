import os

configfile: "sifaka_config.json"

fastq_directory = "/home/thwebste/Data/sifaka"
temp_directory = "temp/"

mmul_path = "reference/GCF_000772875.2_Mmul_8.0.1_genomic.fna"
mmul_prefix = "reference/GCF_000772875.2_Mmul_8.0.1_genomic"
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

all_fastq_prefixes = config["sifaka_fastq_prefixes"] + config["macaque_fastq_prefixes"]
all_samples = config["sifaka_males"] + config["sifaka_females"] + config["macaque_males"] + config["macaque_females"]
sifaka_samples = config["sifaka_males"] + config["sifaka_females"]
macaque_samples = config["macaque_males"] + config["macaque_females"]

rule all:
	input:
		expand(
			"fastqc/{fq_prefix}_fastqc.html", fq_prefix=all_fastq_prefixes),
		expand(
			"adapters/{sample}.adapters.fa", sample=all_samples),
		expand(
			"stats/{sample}.pcoq.{sampling}.mapq.stats",
			sample=sifaka_samples, sampling=["downsampled", "unsampled"]),
		expand(
			"stats/{sample}.hg38.{sampling}.mapq.stats",
			sample=all_samples, sampling=["downsampled", "unsampled"]),
		expand(
			"stats/{sample}.mmul.{sampling}.mapq.stats",
			sample=macaque_samples, sampling=["downsampled", "unsampled"]),
		expand(
			"stats/{sample}.rhemac2.{sampling}.mapq.stats",
			sample=macaque_samples, sampling=["downsampled", "unsampled"]),
		expand(
			"stats/{sample}.pcoq.{sampling}.mapq.stats",
			sample=sifaka_samples, sampling=["downsampled", "unsampled"]),
		expand(
			"stats/{sample}.hg38.{sampling}.mapq.stats",
			sample=all_samples, sampling=["downsampled", "unsampled"]),
		expand(
			"stats/{sample}.mmul.{sampling}.mapq.stats",
			sample=macaque_samples, sampling=["downsampled", "unsampled"]),
		expand(
			"stats/{sample}.rhemac2.{sampling}.mapq.stats",
			sample=macaque_samples, sampling=["downsampled", "unsampled"]),
		expand(
			"callable_sites/combined.{species}.hg38.{chrom}.CHROMcallablesites.{sampling}.bed",
			species=["macaque", "sifaka"], chrom=config["hg38_chroms"],
			sampling=["downsampled", "unsampled"]),
		expand(
			"vcf/sifakas.pcoq.gatk.{sampling}.raw.vcf.gz.tbi",
			sampling=["downsampled", "unsampled"]),
		expand(
			"vcf/macaques.rhemac2.gatk.{sampling}.raw.vcf.gz.tbi",
			sampling=["downsampled", "unsampled"]),
		expand(
			"vcf/macaques.mmul.gatk.{sampling}.raw.vcf.gz.tbi",
			sampling=["downsampled", "unsampled"]),
		expand(
			"vcf/{sample}.macaque.hg38.{chrom}.{sampling}.g.vcf.gz",
			sample=macaque_samples, chrom=config["hg38_chroms"],
			sampling=["downsampled", "unsampled"]),
		expand(
			"vcf/{sample}.sifaka.hg38.{chrom}.{sampling}.g.vcf.gz",
			sample=sifaka_samples, chrom=config["hg38_chroms"],
			sampling=["downsampled", "unsampled"])

		# expand(
		# 	"fastqc/{fq_prefix}_fastqc.html", fq_prefix=all_fastq_prefixes),
		# expand(
		# 	"adapters/{sample}.adapters.fa", sample=all_samples),
		# expand(
		# 	"processed_bams/{sample}.hg38.sorted.mkdup.bam", sample=all_samples),
		# expand(
		# 	"processed_bams/{sample}.pcoq.sorted.mkdup.bam", sample=sifaka_samples),
		# expand(
		# 	"processed_bams/{sample}.mmul.sorted.mkdup.bam", sample=macaque_samples),
		# expand(
		# 	"stats/{sample}.pcoq.sorted.mkdup.bam.stats", sample=sifaka_samples),
		# expand(
		# 	"stats/{sample}.mmul.sorted.mkdup.bam.stats", sample=macaque_samples),
		# expand(
		# 	"stats/{sample}.hg38.sorted.mkdup.bam.stats", sample=all_samples),
		# expand(
		# 	"callable_sites/{sample}.hg38.ONLYcallablesites.bed", sample=all_samples),
		# expand(
		# 	"callable_sites/{sample}.pcoq.ONLYcallablesites.bed", sample=sifaka_samples),
		# expand(
		# 	"callable_sites/{sample}.mmul.ONLYcallablesites.bed", sample=macaque_samples),
		# expand(
		# 	"stats/{sample}.pcoq.mapq.stats", sample=sifaka_samples),
		# expand(
		# 	"stats/{sample}.hg38.mapq.stats", sample=all_samples),
		# expand(
		# 	"stats/{sample}.mmul.mapq.stats", sample=macaque_samples),
		# expand(
		# 	"vcf/sifakas.pcoq.{caller}.raw.vcf.gz.tbi",
		# 	caller=["freebayes", "platypus", "gatk"]),
		# expand(
		# 	"vcf/macaques.mmul.{caller}.raw.vcf.gz.tbi",
		# 	caller=["freebayes", "platypus", "gatk"]),
		# expand(
		# 	"vcf/sifakas.hg38.{caller}.raw.vcf.gz.tbi",
		# 	caller=["freebayes", "platypus", "gatk"]),
		# expand(
		# 	"vcf/macaques.hg38.{caller}.raw.vcf.gz.tbi",
		# 	caller=["freebayes", "platypus", "gatk"])


rule prepare_reference_pcoq_1:
	input:
		config["genome_paths"]["pcoq"]
	output:
		fai = config["genome_paths"]["pcoq"] + ".fai",
		amb = config["genome_paths"]["pcoq"] + ".amb",
		dict = config["genome_prefixes"]["pcoq"] + ".dict"
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
		config["genome_paths"]["hg38"]
	output:
		fai = config["genome_paths"]["hg38"] + ".fai",
		amb = config["genome_paths"]["hg38"] + ".amb",
		dict = config["genome_prefixes"]["hg38"] + ".dict"
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

rule prepare_reference_mmul:
	input:
		config["genome_paths"]["mmul"]
	output:
		fai = config["genome_paths"]["mmul"] + ".fai",
		amb = config["genome_paths"]["mmul"] + ".amb",
		dict = config["genome_prefixes"]["mmul"] + ".dict"
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

rule prepare_reference_rhemac2:
	input:
		config["genome_paths"]["rhemac2"]
	output:
		fai = config["genome_paths"]["rhemac2"] + ".fai",
		amb = config["genome_paths"]["rhemac2"] + ".amb",
		dict = config["genome_prefixes"]["rhemac2"] + ".dict"
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

rule map_and_process_trimmed_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		fai = lambda wildcards: config["genome_paths"][wildcards.genome] + ".fai",
		ref = lambda wildcards: config["genome_paths"][wildcards.genome]
	output:
		"processed_bams/{sample}.{genome}.sorted.mkdup.unsampled.bam"
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

rule downsample_bams:
	input:
		"processed_bams/{sample}.{genome}.sorted.mkdup.unsampled.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.mkdup.downsampled.bam"
	params:
		samtools = samtools_path,
		downsample_fraction = lambda wildcards: config[wildcards.genome][wildcards.sample]
	shell:
		"{params.samtools} view -F 1024 -s 0.{params.downsample_fraction} -b {input} > {output}"

rule index_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam.bai"
	shell:
		"samtools index {input}"

rule bam_stats:
	input:
		"processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam"
	output:
		"stats/{sample}.{genome}.sorted.mkdup.bam.{sampling}.stats"
	shell:
		"samtools stats {input} | grep ^SN | cut -f 2- > {output}"

rule generate_callable_sites:
	input:
		ref = lambda wildcards: config["genome_paths"][wildcards.genome],
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam.bai"
	output:
		"callable_sites/{sample}.{genome}.{sampling}.callablesites"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk,
		summary = "stats/{sample}.{genome}.{sampling}.callable.summary"
	shell:
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T CallableLoci -R {input.ref} -I {input.bam} --minDepth 4 --summary {params.summary} -o {output}"

rule extract_callable_sites:
	input:
		"callable_sites/{sample}.{genome}.{sampling}.callablesites"
	output:
		"callable_sites/{sample}.{genome}.ONLYcallablesites.{sampling}.bed"
	shell:
		"sed -e '/CALLABLE/!d' {input} > {output}"

rule combine_callable_sites_rhemac2:
	input:
		expand(
			"callable_sites/{sample}.rhemac2.ONLYcallablesites.{{sampling}}.bed",
			sample=macaque_samples)
	output:
		"callable_sites/combined.rhemac2.COMBINEDcallablesites.{sampling}.bed"
	params:
		bedtools = bedtools_path
	shell:
		"cat {input} | sort -k1,1 -k2,2n | {params.bedtools} merge -i stdin > {output}"

rule combine_callable_sites_pcoq:
	input:
		expand(
			"callable_sites/{sample}.pcoq.ONLYcallablesites.{{sampling}}.bed",
			sample=sifaka_samples)
	output:
		"callable_sites/combined.pcoq.COMBINEDcallablesites.{sampling}.bed"
	params:
		bedtools = bedtools_path
	shell:
		"cat {input} | sort -k1,1 -k2,2n | {params.bedtools} merge -i stdin > {output}"

rule combine_callable_sites_hg38_sifakas:
	input:
		expand(
			"callable_sites/{sample}.hg38.ONLYcallablesites.{{sampling}}.bed",
			sample=sifaka_samples)
	output:
		"callable_sites/combined.sifaka.hg38.COMBINEDcallablesites.{sampling}.bed"
	params:
		bedtools = bedtools_path
	shell:
		"cat {input} | sort -k1,1 -k2,2n | {params.bedtools} merge -i stdin > {output}"

rule combine_callable_sites_hg38_macaques:
	input:
		expand(
			"callable_sites/{sample}.hg38.ONLYcallablesites.{{sampling}}.bed",
			sample=macaque_samples)
	output:
		"callable_sites/combined.macaque.hg38.COMBINEDcallablesites.{sampling}.bed"
	params:
		bedtools = bedtools_path
	shell:
		"cat {input} | sort -k1,1 -k2,2n | {params.bedtools} merge -i stdin > {output}"

rule split_hg38_callable_bed_by_chrom:
	input:
		callable = "callable_sites/combined.{species}.hg38.COMBINEDcallablesites.{sampling}.bed"
	output:
		"callable_sites/combined.{species}.hg38.{chrom}.CHROMcallablesites.{sampling}.bed"
	shell:
		"""awk '$1=="{wildcards.chrom}"' {input.callable} > {output}"""

rule mapq_check:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam.bai"
	output:
		"stats/{sample}.{genome}.{sampling}.mapq.stats"
	threads: 4
	shell:
		"scripts/mapqs -infile {input.bam} -outfile {output} -threads {threads}"

rule gatk_gvcf_pcoq:
	input:
		ref = config["genome_paths"]["pcoq"],
		bam = "processed_bams/{sample}.pcoq.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.pcoq.sorted.mkdup.{sampling}.bam.bai",
		callable = "callable_sites/combined.pcoq.COMBINEDcallablesites.{sampling}.bed"
	output:
		"vcf/{sample}.pcoq.{sampling}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads: 4
	shell:
		"java -Xmx16g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.callable} --emitRefConfidence GVCF -o {output}"

rule gatk_gvcf_rhemac2:
	input:
		ref = config["genome_paths"]["rhemac2"],
		bam = "processed_bams/{sample}.rhemac2.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.rhemac2.sorted.mkdup.{sampling}.bam.bai",
		callable = "callable_sites/combined.rhemac2.COMBINEDcallablesites.{sampling}.bed"
	output:
		"vcf/{sample}.rhemac2.{sampling}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads: 4
	shell:
		"java -Xmx16g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.callable} --emitRefConfidence GVCF -o {output}"

rule gatk_gvcf_mmul:
	input:
		ref = config["genome_paths"]["mmul"],
		bam = "processed_bams/{sample}.mmul.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.mmul.sorted.mkdup.{sampling}.bam.bai",
		callable = "callable_sites/combined.mmul.COMBINEDcallablesites.{sampling}.bed"
	output:
		"vcf/{sample}.mmul.{sampling}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads: 4
	shell:
		"java -Xmx16g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.callable} --emitRefConfidence GVCF -o {output}"

#########################

rule gatk_gvcf_hg38_sifaka:
	input:
		ref = config["genome_paths"]["hg38"],
		bam = "processed_bams/{sample}.hg38.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.hg38.sorted.mkdup.{sampling}.bam.bai",
		callable = "callable_sites/combined.sifaka.hg38.ONLYcallablesites.bed"
	output:
		"vcf/{sample}.sifaka.hg38.{chrom}.{sampling}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk,
		chromosome = "{chrom}"
	threads: 4
	shell:
		"java -Xmx16g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.callable} -L {params.chromosome} --emitRefConfidence GVCF -o {output}"

rule gatk_gvcf_hg38_macaque:
	input:
		ref = config["genome_paths"]["hg38"],
		bam = "processed_bams/{sample}.hg38.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.hg38.sorted.mkdup.{sampling}.bam.bai",
		callable = "callable_sites/combined.macaque.hg38.ONLYcallablesites.bed"
	output:
		"vcf/{sample}.macaque.hg38.{chrom}.{sampling}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk,
		chromosome = "{chrom}"
	threads: 4
	shell:
		"java -Xmx16g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.callable} -L {params.chromosome} --emitRefConfidence GVCF -o {output}"

rule gatk_cat_variants_hg38:
	input:
		ref = hg38_path,
		gvcfs = expand(
			"vcf/{{sample}}.{{species}}.hg38.{chrom}.g.vcf.gz", chrom=config["hg38_chroms"])
	output:
		"vcf/{sample}.{species}.hg38.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("-V " + i)
		variant_files = " ".join(variant_files)
		shell("java -Xmx16g -Djava.io.tmpdir={params.temp_dir} -cp {params.gatk_path} org.broadinstitute.gatk.tools.CatVariants -R {input.ref} -o {output}")

rule genotype_gvcfs_pcoq:
	input:
		ref = pcoq_1_path,
		gvcfs = expand("vcf/{sample}.pcoq.g.vcf.gz", sample=sifaka_samples)
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

rule genotype_gvcfs_macaque:
	input:
		ref = mmul_path,
		gvcfs = expand("vcf/{sample}.mmul.g.vcf.gz", sample=macaque_samples)
	output:
		v = "vcf/macaques.mmul.gatk.raw.vcf.gz"
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

rule genotype_gvcfs_hg38_sifaka:
	input:
		ref = hg38_path,
		gvcfs = expand("vcf/{sample}.sifaka.hg38.g.vcf.gz", sample=sifaka_samples)
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

rule genotype_gvcfs_hg38_macaque:
	input:
		ref = hg38_path,
		gvcfs = expand("vcf/{sample}.macaque.hg38.g.vcf.gz", sample=macaque_samples)
	output:
		v = "vcf/macaques.hg38.gatk.raw.vcf.gz"
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

rule zip_vcf:
	input:
		vcf = "vcf/{species}.{genome}.{caller}.raw.vcf"
	output:
		"vcf/{species}.{genome}.{caller}.raw.vcf.gz"
	params:
		bgzip = bgzip_path
	shell:
		"{params.bgzip} {input.vcf}"

rule index_zipped_vcf:
	input:
		vcf = "vcf/{species}.{genome}.{caller}.raw.vcf.gz"
	output:
		"vcf/{species}.{genome}.{caller}.raw.vcf.gz.tbi"
	params:
		tabix = tabix_path
	shell:
		"{params.tabix} -p vcf {input}"
