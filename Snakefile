import os

configfile: "sifaka_config.json"

fastq_directory = "/home/thwebste/Data/sifaka"
temp_directory = "temp/"

bbduksh_path = "bbduk.sh"
bbmerge_sh_path = "bbmerge.sh"
bcftools_path = "bcftools"
bedops_path = "bedops"
bedtools_path = "bedtools"
bgzip_path = "bgzip"
bwa_path = "bwa"
fastqc_path = "fastqc"
freebayes_path = "freebayes"
gatk = "/home/thwebste/Tools/GenomeAnalysisTK_37.jar"
gff2bed_path = "gff2bed"
multiqc_path = "multiqc"
platypus_path = "platypus"
samblaster_path = "samblaster"
samtools_path = "samtools"
sort_bed_path = "sort-bed"
tabix_path = "tabix"
vep_path = "vep"
# xyalign_path =

all_fastq_prefixes = config["sifaka_fastq_prefixes"] + config["macaque_fastq_prefixes"]
all_samples_main = config["sifaka_males"] + config["sifaka_females"] + config["macaque_males"] + config["macaque_females"]
sifaka_samples = config["sifaka_males"] + config["sifaka_females"]
macaque_samples = config["macaque_males"] + config["macaque_females"]

# Added second set of sifakas for capture kit comparison
secondary_sifaka_samples = config["secondary_sifakas"]
combined_sifaka_samples = sifaka_samples + secondary_sifaka_samples
all_samples_both_sifakas = all_samples_main + secondary_sifaka_samples

rule all:
	input:
		"multiqc/multiqc_report.html",
		expand(
			"adapters/{sample}.adapters.fa", sample=all_samples_both_sifakas),
		"multiqc_trimmed/multiqc_report.html",
		expand(
			"stats/{sample}.pcoq.{sampling}.mapq.stats",
			sample=combined_sifaka_samples,
			sampling=["downsampled", "unsampled"]),
		expand(
			"stats/{sample}.hg38.{sampling}.mapq.stats",
			sample=all_samples_both_sifakas,
			sampling=["downsampled", "unsampled"]),
		expand(
			"stats/{sample}.mmul.{sampling}.mapq.stats",
			sample=macaque_samples,
			sampling=["downsampled", "unsampled"]),
		expand(
			"stats/{sample}.pcoq.sorted.mkdup.bam.{sampling}.{suffix}",
			sample=combined_sifaka_samples,
			sampling=["downsampled", "unsampled"],
			suffix=["stats", "nodup.properpair.stats"]),
		expand(
			"stats/{sample}.hg38.sorted.mkdup.bam.{sampling}.{suffix}",
			sample=all_samples_both_sifakas,
			sampling=["downsampled", "unsampled"],
			suffix=["stats", "nodup.properpair.stats"]),
		expand(
			"stats/{sample}.mmul.sorted.mkdup.bam.{sampling}.{suffix}",
			sample=macaque_samples,
			sampling=["downsampled", "unsampled"],
			suffix=["stats", "nodup.properpair.stats"]),
		expand(
			"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.cdscoverage.bed",
			sample=combined_sifaka_samples, genome=["hg38", "pcoq"],
			sampling=["downsampled", "unsampled"]),
		expand(
			"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.cdscoverage.bed",
			sample=macaque_samples, genome=["hg38", "mmul"],
			sampling=["downsampled", "unsampled"]),
		expand(
			"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.{region}.hist",
			sample=combined_sifaka_samples, genome=["hg38", "pcoq"],
			sampling=["downsampled", "unsampled"],
			region=["cds", "exon", "gene", "utr", "intron", "intergenic"]),
		expand(
			"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.{region}.hist",
			sample=macaque_samples, genome=["hg38", "mmul"],
			sampling=["downsampled", "unsampled"],
			region=["cds", "exon", "gene", "utr", "intron", "intergenic"]),
		expand(
			"vcf/sifakas.pcoq.{caller}.{sampling}.filtered.vep.vcf.gz",
			caller=["gatk", "freebayes"],
			sampling=["downsampled"]),
		expand(
			"vcf/macaques.mmul.{caller}.{sampling}.filtered.vep.vcf.gz",
			caller=["gatk", "freebayes"],
			sampling=["downsampled"]),
		expand(
			"vcf/sifakas.hg38.{caller}.{sampling}.filtered.vep.vcf.gz",
			caller=["gatk", "freebayes"],
			sampling=["downsampled"]),
		expand(
			"vcf/macaques.hg38.{caller}.{sampling}.filtered.vep.vcf.gz",
			caller=["gatk", "freebayes"],
			sampling=["downsampled"]),
		expand(
			"stats/sifakas.{genome}.{caller}.{sampling}.{region}.INTERSECTION.vcf.stats",
			sample=combined_sifaka_samples, genome=["hg38", "pcoq"],
			caller=["gatk", "freebayes"],
			sampling=["downsampled"],
			region=["cds", "exon", "gene", "utr", "intron"]),
		expand(
			"stats/macaques.{genome}.{caller}.{sampling}.{region}.INTERSECTION.vcf.stats",
			sample=macaque_samples, genome=["hg38", "mmul"],
			caller=["gatk", "freebayes"],
			sampling=["downsampled"],
			region=["cds", "exon", "gene", "utr", "intron"]),
		expand(
			"stats/sifakas.{genome}.{caller}.{sampling}.intergenic.SUBTRACTION.vcf.stats",
			sample=combined_sifaka_samples, genome=["hg38", "pcoq"],
			caller=["gatk", "freebayes"],
			sampling=["downsampled"]),
		expand(
			"stats/macaques.{genome}.{caller}.{sampling}.intergenic.SUBTRACTION.vcf.stats",
			sample=macaque_samples, genome=["hg38", "mmul"],
			caller=["gatk", "freebayes"],
			sampling=["downsampled"]),
		"stats/compiled.bcftools_stats.csv"



rule prepare_reference:
	input:
		ref = lambda wildcards: config["genome_paths"][wildcards.assembly]
	output:
		new = "new_reference/{assembly}.fasta",
		fai = "new_reference/{assembly}.fasta.fai",
		amb = "new_reference/{assembly}.fasta.amb",
		dict = "new_reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		shell(
			"ln -s ../{} {{output.new}} && touch -h {{output.new}}".format(input.ref))
		# faidx
		shell(
			"{params.samtools} faidx {output.new}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {output.new}")
		# bwa
		shell(
			"{params.bwa} index {output.new}")

rule get_annotation:
	output:
		"new_reference/{genome}.gff"
	params:
		web_address = lambda wildcards: config["annotation_address"][wildcards.genome],
		initial_output = "new_reference/{genome}.gff.gz"
	run:
		shell("wget {params.web_address} -O {params.initial_output}")
		shell("gunzip {params.initial_output}")

rule extract_exons_from_gff:
	input:
		"new_reference/{genome}.gff"
	output:
		"regions/{genome}.exon.gff"
	shell:
		"""awk '($3 == "exon")' {input} > {output}"""

rule extract_genes_from_gff:
	input:
		"new_reference/{genome}.gff"
	output:
		"regions/{genome}.gene.gff"
	shell:
		"""awk '($3 == "gene")' {input} > {output}"""

rule extract_cds_from_gff:
	input:
		"new_reference/{genome}.gff"
	output:
		"regions/{genome}.cds.gff"
	shell:
		"""awk '($3 == "CDS")' {input} > {output}"""

rule create_introns:
	input:
		exon = "regions/{genome}.exon.gff",
		gene = "regions/{genome}.gene.gff"
	output:
		"regions/{genome}.intron.gff"
	shell:
		"bedtools subtract -a {input.gene} -b {input.exon} > {output}"

rule create_utr:
	input:
		exon = "regions/{genome}.exon.gff",
		cds = "regions/{genome}.cds.gff"
	output:
		"regions/{genome}.utr.gff"
	shell:
		"bedtools subtract -a {input.exon} -b {input.cds} > {output}"

rule convert_gff_coordinates:
	input:
		"regions/{genome}.{region}.gff"
	output:
		"regions/{genome}.{region}.converted.gff"
	run:
		if wildcards.genome == "hg38":
			shell(
				"python scripts/Convert_hg38_names.py --gff {input} --outfile {output}")
		else:
			shell(
				"ln -s ../{} {{output}} && touch -h {{output}}".format(input))

rule fastqc_analysis:
	input:
		os.path.join(fastq_directory, "{fq_prefix}.fastq.gz")
	output:
		"fastqc/{fq_prefix}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc {input}"

rule multiqc_analysis:
	input:
		expand("fastqc/{fq_prefix}_fastqc.html", fq_prefix=all_fastq_prefixes)
	output:
		"multiqc/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} -o multiqc fastqc"

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
		"{params.bbduksh} -Xmx3g in1={input.fq1} in2={input.fq2} out1={output.out_fq1} out2={output.out_fq2} ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe qtrim=rl trimq=10"

rule fastqc_analysis_trimmed_paired:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		html1 = "fastqc_trimmed/{sample}_trimmed_read1_fastqc.html",
		html2 = "fastqc_trimmed/{sample}_trimmed_read2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_trimmed {input.fq1} {input.fq2}"

rule trimmed_multiqc_analysis:
	input:
		expand(
			"fastqc_trimmed/{sample}_trimmed_{read}_fastqc.html",
			sample=all_samples_both_sifakas, read=["read1", "read2"])
	output:
		"multiqc_trimmed/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} -o multiqc_trimmed fastqc_trimmed"

rule map_and_process_trimmed_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		fai = "new_reference/{genome}.fasta.fai",
		ref = "new_reference/{genome}.fasta"
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
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule bam_stats:
	input:
		"processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam"
	output:
		"stats/{sample}.{genome}.sorted.mkdup.bam.{sampling}.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} stats {input} | grep ^SN | cut -f 2- > {output}"

rule bam_stats_on_passing_reads:
	input:
		"processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam"
	output:
		"stats/{sample}.{genome}.sorted.mkdup.bam.{sampling}.nodup.properpair.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} view -F 1024 -f 2 -b {input} | {params.samtools} stats - | grep ^SN | cut -f 2- > {output}"

rule mapq_check:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam.bai"
	output:
		"stats/{sample}.{genome}.{sampling}.mapq.stats"
	threads:
		4
	shell:
		"scripts/mapqs -infile {input.bam} -outfile {output} -threads {threads}"

rule make_bedtools_genome_file:
	input:
		fai = "new_reference/{genome}.fasta.fai"
	output:
		"new_reference/{genome}.genome"
	shell:
		"awk -v OFS='\t' {{'print $1,$2'}} {input.fai} > {output}"

rule genome_cov:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam",
		idx = "processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam.bai",
		genome = "new_reference/{genome}.genome"
	output:
		"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov"
	params:
		samtools = samtools_path,
		bedtools = bedtools_path
	shell:
		"{params.samtools} view {input.bam} -b -F 1024 -q 20 | "
		"{params.bedtools} genomecov -bg -ibam stdin -g {input.genome} > {output}"

rule bedops_sort_genome_cov:
	input:
		"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov"
	output:
		"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.bedopssorted.bed"
	params:
		sort_bed = sort_bed_path
	shell:
		"{params.sort_bed} {input} > {output}"

rule gff2bed:
	input:
		"regions/{genome}.{region}.converted.gff"
	output:
		"regions/{genome}.{region}.converted.bedopssorted.bed"
	params:
		gff2bed = gff2bed_path,
		bedtools = bedtools_path
	shell:
		"cat {input} | {params.gff2bed} | {params.bedtools} merge > {output}"

rule cds_target_coverage:
	input:
		sample = "results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.bedopssorted.bed",
		region = "regions/{genome}.cds.converted.bedopssorted.bed"
	output:
		"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.cdscoverage.bed"
	params:
		bedtools = bedtools_path
	shell:
		"{params.bedtools} coverage -sorted -a {input.region} -b {input.sample} > {output}"

rule bedtools_intersect_regions:
	input:
		sample = "results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.bedopssorted.bed",
		region = "regions/{genome}.{region}.converted.bedopssorted.bed"
	output:
		"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.bedopssorted.{region}.INTERSECTION.bed"
	params:
		bedtools = bedtools_path
	shell:
		"{params.bedtools} intersect -a {input.sample} -b {input.region} > {output}"

rule bedtools_find_intergenic:
	input:
		sample = "results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.bedopssorted.bed",
		region = "regions/{genome}.gene.converted.bedopssorted.bed"
	output:
		"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.bedopssorted.intergenic.SUBTRACTION.bed"
	params:
		bedtools = bedtools_path
	shell:
		"{params.bedtools} subtract -a {input.sample} -b {input.region} > {output}"

rule compute_histogram_from_bed:
	input:
		"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.bedopssorted.{region}.INTERSECTION.bed"
	output:
		"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.{region}.hist"
	shell:
		"python scripts/Compute_histogram_from_bed.py --bed {input} --outfile {output}"

rule compute_histogram_from_bed_intergenic:
	input:
		"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.bedopssorted.intergenic.SUBTRACTION.bed"
	output:
		"results/{sample}.{genome}.{sampling}.mapq20_noDup.genome_cov.intergenic.hist"
	shell:
		"python scripts/Compute_histogram_from_bed.py --bed {input} --outfile {output}"

rule generate_callable_sites:
	input:
		ref = "new_reference/{genome}.fasta",
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

rule combine_callable_sites_mmul:
	input:
		expand(
			"callable_sites/{sample}.mmul.ONLYcallablesites.{{sampling}}.bed",
			sample=macaque_samples)
	output:
		"callable_sites/combined.mmul.COMBINEDcallablesites.{sampling}.bed"
	params:
		bedtools = bedtools_path
	shell:
		"cat {input} | sort -k1,1 -k2,2n | {params.bedtools} merge -i stdin > {output}"

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

rule gatk_gvcf_hg38_sifaka:
	input:
		ref = "new_reference/hg38.fasta",
		bam = "processed_bams/{sample}.hg38.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.hg38.sorted.mkdup.{sampling}.bam.bai",
		callable = "callable_sites/combined.sifaka.hg38.COMBINEDcallablesites.{sampling}.bed"
	output:
		"vcf/{sample}.sifaka.hg38.{sampling}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads:
		8
	shell:
		"java -Xmx30g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.callable} --emitRefConfidence GVCF -o {output}"

rule gatk_gvcf_hg38_macaques:
	input:
		ref = "new_reference/hg38.fasta",
		bam = "processed_bams/{sample}.hg38.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.hg38.sorted.mkdup.{sampling}.bam.bai",
		callable = "callable_sites/combined.macaque.hg38.COMBINEDcallablesites.{sampling}.bed"
	output:
		"vcf/{sample}.macaque.hg38.{sampling}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads:
		8
	shell:
		"java -Xmx30g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.callable} --emitRefConfidence GVCF -o {output}"

rule gatk_gvcf_pcoq:
	input:
		ref = "new_reference/pcoq.fasta",
		bam = "processed_bams/{sample}.pcoq.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.pcoq.sorted.mkdup.{sampling}.bam.bai",
		callable = "callable_sites/combined.pcoq.COMBINEDcallablesites.{sampling}.bed"
	output:
		"vcf/{sample}.pcoq.{sampling}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads:
		8
	shell:
		"java -Xmx30g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.callable} --emitRefConfidence GVCF -o {output}"

rule gatk_gvcf_mmul:
	input:
		ref = "new_reference/mmul.fasta",
		bam = "processed_bams/{sample}.mmul.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.mmul.sorted.mkdup.{sampling}.bam.bai",
		callable = "callable_sites/combined.mmul.COMBINEDcallablesites.{sampling}.bed"
	output:
		"vcf/{sample}.mmul.{sampling}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads: 8
	shell:
		"java -Xmx30g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.callable} --emitRefConfidence GVCF -o {output}"

rule genotype_gvcfs_hg38_sifaka:
	input:
		ref = "new_reference/hg38.fasta",
		gvcfs = expand(
			"vcf/{sample}.sifaka.hg38.{{sampling}}.g.vcf.gz",
			sample=sifaka_samples)
	output:
		v = "vcf/macaques.hg38.gatk.{sampling}.raw.vcf"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell("java -Xmx16g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T GenotypeGVCFs -R {input.ref} {variant_files} -o {output.v} --includeNonVariantSites")

rule genotype_gvcfs_hg38_macaque:
	input:
		ref = "new_reference/hg38.fasta",
		gvcfs = expand(
			"vcf/{sample}.macaque.hg38.{{sampling}}.g.vcf.gz",
			sample=macaque_samples)
	output:
		v = "vcf/sifakas.hg38.gatk.{sampling}.raw.vcf"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads:
		8
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell("java -Xmx30g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T GenotypeGVCFs -R {input.ref} {variant_files} -o {output.v} --includeNonVariantSites")

rule genotype_gvcfs_pcoq:
	input:
		ref = "new_reference/pcoq.fasta",
		gvcfs = expand(
			"vcf/{sample}.pcoq.{{sampling}}.g.vcf.gz",
			sample=sifaka_samples)
	output:
		v = "vcf/sifakas.pcoq.gatk.{sampling}.raw.vcf"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads:
		8
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell("java -Xmx30g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T GenotypeGVCFs -R {input.ref} {variant_files} -o {output.v} --includeNonVariantSites")

rule genotype_gvcfs_mmul:
	input:
		ref = "new_reference/mmul.fasta",
		gvcfs = expand(
			"vcf/{sample}.mmul.{{sampling}}.g.vcf.gz",
			sample=macaque_samples)
	output:
		v = "vcf/macaques.mmul.gatk.{sampling}.raw.vcf"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk
	threads:
		8
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell("java -Xmx30g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T GenotypeGVCFs -R {input.ref} {variant_files} -o {output.v} --includeNonVariantSites")

rule freebayes_call_mmul:
	input:
		ref = "new_reference/mmul.fasta",
		bams = expand(
			"processed_bams/{sample}.mmul.sorted.mkdup.{{sampling}}.bam",
			sample=macaque_samples),
		bais = expand(
			"processed_bams/{sample}.mmul.sorted.mkdup.{{sampling}}.bam.bai",
			sample=macaque_samples),
		callable = "callable_sites/combined.mmul.COMBINEDcallablesites.{sampling}.bed"
	output:
		vcf = "vcf/macaques.mmul.freebayes.{sampling}.raw.vcf"
	params:
		freebayes = freebayes_path
	shell:
		"{params.freebayes} -f {input.ref} -v {output} --genotype-qualities --targets {input.callable} {input.bams}"

rule freebayes_call_pcoq:
	input:
		ref = "new_reference/pcoq.fasta",
		bams = expand(
			"processed_bams/{sample}.pcoq.sorted.mkdup.{{sampling}}.bam",
			sample=sifaka_samples),
		bais = expand(
			"processed_bams/{sample}.pcoq.sorted.mkdup.{{sampling}}.bam.bai",
			sample=sifaka_samples),
		callable = "callable_sites/combined.pcoq.COMBINEDcallablesites.{sampling}.bed"
	output:
		vcf = "vcf/sifakas.pcoq.freebayes.{sampling}.raw.vcf"
	params:
		freebayes = freebayes_path
	shell:
		"{params.freebayes} -f {input.ref} -v {output} --genotype-qualities --targets {input.callable} {input.bams}"

rule freebayes_call_hg38_sifakas:
	input:
		ref = "new_reference/hg38.fasta",
		bams = expand(
			"processed_bams/{sample}.hg38.sorted.mkdup.{{sampling}}.bam",
			sample=sifaka_samples),
		bais = expand(
			"processed_bams/{sample}.hg38.sorted.mkdup.{{sampling}}.bam.bai",
			sample=sifaka_samples),
		callable = "callable_sites/combined.sifaka.hg38.COMBINEDcallablesites.{sampling}.bed"
	output:
		vcf = "vcf/sifakas.hg38.freebayes.{sampling}.raw.vcf"
	params:
		freebayes = freebayes_path
	shell:
		"{params.freebayes} -f {input.ref} -v {output} --genotype-qualities --targets {input.callable} {input.bams}"

rule freebayes_call_hg38_macaques:
	input:
		ref = "new_reference/hg38.fasta",
		bams = expand(
			"processed_bams/{sample}.hg38.sorted.mkdup.{{sampling}}.bam",
			sample=macaque_samples),
		bais = expand(
			"processed_bams/{sample}.hg38.sorted.mkdup.{{sampling}}.bam.bai",
			sample=macaque_samples),
		callable = "callable_sites/combined.macaque.hg38.COMBINEDcallablesites.{sampling}.bed"
	output:
		vcf = "vcf/macaques.hg38.freebayes.{sampling}.raw.vcf"
	params:
		freebayes = freebayes_path
	shell:
		"{params.freebayes} -f {input.ref} -v {output} --genotype-qualities --targets {input.callable} {input.bams}"

rule zip_vcf:
	input:
		vcf = "vcf/{species}.{genome}.{caller}.{sampling}.raw.vcf"
	output:
		"vcf/{species}.{genome}.{caller}.{sampling}.raw.vcf.gz"
	params:
		bgzip = bgzip_path
	shell:
		"{params.bgzip} {input.vcf}"

rule index_zipped_vcf:
	input:
		vcf = "vcf/{species}.{genome}.{caller}.{sampling}.raw.vcf.gz"
	output:
		"vcf/{species}.{genome}.{caller}.{sampling}.raw.vcf.gz.tbi"
	params:
		tabix = tabix_path
	shell:
		"{params.tabix} -p vcf {input}"

rule filter_vcfs:
	input:
		vcf = "vcf/{species}.{genome}.{caller}.{sampling}.raw.vcf.gz",
		tbi = "vcf/{species}.{genome}.{caller}.{sampling}.raw.vcf.gz.tbi"
	output:
		"vcf/{species}.{genome}.{caller}.{sampling}.filtered.vcf"
	params:
		filter_script = "scripts/Filter_vcf.py",
		num_samples = 4,
		variant_caller = "{caller}"
	shell:
		"python {params.filter_script} --vcf {input.vcf} --output_vcf {output} --variant_caller {params.variant_caller} --QUAL 30 --sample_depth 8 --min_samples {params.num_samples} --min_support 3 --genotype_quality 30"

rule zip_filtered_vcfs:
	input:
		vcf = "vcf/{species}.{genome}.{caller}.{sampling}.filtered.vcf"
	output:
		"vcf/{species}.{genome}.{caller}.{sampling}.filtered.vcf.gz"
	shell:
		"bgzip {input.vcf}"

rule index_zipped_filtered_vcfs:
	input:
		"vcf/{species}.{genome}.{caller}.{sampling}.filtered.vcf.gz"
	output:
		"vcf/{species}.{genome}.{caller}.{sampling}.filtered.vcf.gz.tbi"
	shell:
		"tabix -p vcf {input}"

rule download_hg38_cache:
	output:
		"new_reference/homo_sapiens"
	params:
		web_address = lambda wildcards: config["vep_address"]["hg38"],
		initial_output = "new_reference/homo_sapiens_vep_90_GRCh38.tar.gz"
	run:
		shell("wget {params.web_address} -O {params.initial_output}")
		shell("tar xzf {params.initial_output} -C new_reference")

rule bgzip_gff_files:
	input:
		"new_reference/{genome}.gff"
	output:
		"new_reference/{genome}.gff.gz"
	params:
		bgzip = bgzip_path
	shell:
		"""grep -v "#" {input} | """
		"""sort -k1,1 -k4,4n -k5,5n -t$'\t' | {params.bgzip} -c > {output}"""

rule index_bgzipped_gff_files:
	input:
		"new_reference/{genome}.gff.gz"
	output:
		"new_reference/{genome}.gff.gz.tbi"
	params:
		tabix = tabix_path
	shell:
		"{params.tabix} -p gff {input}"

rule vep_annotation:
	input:
		ref = "new_reference/{genome}.fasta",
		gff = "new_reference/{genome}.gff.gz",
		gff_idx = "new_reference/{genome}.gff.gz.tbi",
		vcf = "vcf/{species}.{genome}.{caller}.{sampling}.filtered.vcf.gz",
		idx = "vcf/{species}.{genome}.{caller}.{sampling}.filtered.vcf.gz.tbi",
		hg38_cache = "new_reference/homo_sapiens"
	output:
		vcf = "vcf/{species}.{genome}.{caller}.{sampling}.filtered.vep.vcf.gz",
		html = "vcf/{species}.{genome}.{caller}.{sampling}.filtered.vep.vcf.gz_summary.html"
	params:
		vep = vep_path,
		genome = "{genome}",
		cache_base = "new_reference/"
	run:
		if params.genome == "hg38":
			shell(
				"{params.vep} -i {input.vcf} --dir_cache {params.cache_base} "
				"-o {output.vcf} --compress_output bgzip --cache --offline --cache_version 90 --vcf")
		else:
			shell(
				"{params.vep} -i {input.vcf} -gff {input.gff} -fasta {input.ref} "
				"-o {output.vcf} --compress_output bgzip --vcf")

rule vcf_intersect_by_region:
	input:
		vcf = "vcf/{species}.{genome}.{caller}.{sampling}.filtered.vcf.gz",
		region = "regions/{genome}.{region}.converted.bedopssorted.bed"
	output:
		"vcf_regions/{species}.{genome}.{caller}.{sampling}.{region}.INTERSECTION.filtered.vcf.gz"
	params:
		bedtools = bedtools_path,
		bgzip = bgzip_path
	shell:
		"{params.bedtools} intersect -header -a {input.vcf} -b {input.region} | "
		"{params.bgzip} > {output}"

rule vcf_intersect_intergenic:
	input:
		vcf = "vcf/{species}.{genome}.{caller}.{sampling}.filtered.vcf.gz",
		region = "regions/{genome}.gene.converted.bedopssorted.bed"
	output:
		"vcf_regions/{species}.{genome}.{caller}.{sampling}.intergenic.SUBTRACTION.filtered.vcf.gz"
	params:
		bedtools = bedtools_path,
		bgzip = bgzip_path
	shell:
		"{params.bedtools} subtract -header -a {input.vcf} -b {input.region} | "
		"{params.bgzip} > {output}"

rule index_zipped_vcf_by_region:
	input:
		vcf = "vcf_regions/{species}.{genome}.{caller}.{sampling}.{region}.INTERSECTION.filtered.vcf.gz"
	output:
		"vcf_regions/{species}.{genome}.{caller}.{sampling}.{region}.INTERSECTION.filtered.vcf.gz.tbi"
	params:
		tabix = tabix_path
	shell:
		"{params.tabix} -p vcf {input}"

rule index_zipped_vcf_intergenic:
	input:
		vcf = "vcf_regions/{species}.{genome}.{caller}.{sampling}.intergenic.SUBTRACTION.filtered.vcf.gz"
	output:
		"vcf_regions/{species}.{genome}.{caller}.{sampling}.intergenic.SUBTRACTION.filtered.vcf.gz.tbi"
	params:
		tabix = tabix_path
	shell:
		"{params.tabix} -p vcf {input}"

rule vcf_stats_by_region:
	input:
		vcf = "vcf_regions/{species}.{genome}.{caller}.{sampling}.{region}.INTERSECTION.filtered.vcf.gz",
		tbi = "vcf_regions/{species}.{genome}.{caller}.{sampling}.{region}.INTERSECTION.filtered.vcf.gz.tbi"
	output:
		"stats/{species}.{genome}.{caller}.{sampling}.{region}.INTERSECTION.vcf.stats"
	params:
		bcftools = bcftools_path
	shell:
		"{params.bcftools} stats {input.vcf} | grep ^SN > {output}"

rule vcf_stats_intergenic:
	input:
		vcf = "vcf_regions/{species}.{genome}.{caller}.{sampling}.intergenic.SUBTRACTION.filtered.vcf.gz",
		tbi = "vcf_regions/{species}.{genome}.{caller}.{sampling}.intergenic.SUBTRACTION.filtered.vcf.gz.tbi"
	output:
		"stats/{species}.{genome}.{caller}.{sampling}.intergenic.SUBTRACTION.vcf.stats"
	params:
		bcftools = bcftools_path
	shell:
		"{params.bcftools} stats {input.vcf} | grep ^SN > {output}"

rule compile_vcf_stats:
	input:
		s_region = expand(
			"stats/sifakas.{genome}.{caller}.{sampling}.{region}.INTERSECTION.vcf.stats",
			sample=combined_sifaka_samples, genome=["hg38", "pcoq"],
			caller=["gatk", "freebayes"],
			sampling=["downsampled"],
			region=["cds", "exon", "gene", "utr", "intron"]),
		m_region = expand(
			"stats/macaques.{genome}.{caller}.{sampling}.{region}.INTERSECTION.vcf.stats",
			sample=macaque_samples, genome=["hg38", "mmul"],
			caller=["gatk", "freebayes"],
			sampling=["downsampled"],
			region=["cds", "exon", "gene", "utr", "intron"]),
		s_inter = expand(
			"stats/sifakas.{genome}.{caller}.{sampling}.intergenic.SUBTRACTION.vcf.stats",
			sample=combined_sifaka_samples, genome=["hg38", "pcoq"],
			caller=["gatk", "freebayes"],
			sampling=["downsampled"]),
		m_inter = expand(
			"stats/macaques.{genome}.{caller}.{sampling}.intergenic.SUBTRACTION.vcf.stats",
			sample=macaque_samples, genome=["hg38", "mmul"],
			caller=["gatk", "freebayes"],
			sampling=["downsampled"])
	output:
		"stats/compiled.bcftools_stats.csv"
	shell:
		"python scripts/Compile_bcftools_stats_output.py --input_files {input} "
		"--output_file {output}"

rule create_coverage_histograms:
	input:
		regions = "regions/{genome}.{region}.converted.gff",
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam.bai"
	output:
		"coverage/{sample}.{genome}.{sampling}.{region}.hist.txt"
	params:
		bedtools = bedtools_path
	shell:
		"{params.bedtools} coverage -hist -abam {input.bam} -b {input.regions} | grep ^all > {output}"

rule generate_callable_sites_analysis:
	input:
		ref = "new_reference/{genome}.fasta",
		bam = "processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.mkdup.{sampling}.bam.bai"
	output:
		"callable_sites/{sample}.{genome}.{sampling}.callablesitesFORANALYSIS"
	params:
		temp_dir = temp_directory,
		gatk_path = gatk,
		summary = "stats/{sample}.{genome}.{sampling}.callableFORANALYSIS.summary"
	shell:
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T CallableLoci -R {input.ref} -I {input.bam} --minDepth 8 --summary {params.summary} -o {output}"

rule extract_callable_sites_analysis:
	input:
		"callable_sites/{sample}.{genome}.{sampling}.callablesitesFORANALYSIS"
	output:
		"callable_sites/{sample}.{genome}.{sampling}.ONLYcallablesitesFORANALYSIS.bed"
	shell:
		"sed -e '/CALLABLE/!d' {input} > {output}"

rule bedops_sort_callable:
	input:
		"callable_sites/{sample}.{genome}.{sampling}.ONLYcallablesitesFORANALYSIS.bed"
	output:
		"callable_sites/{sample}.{genome}.ONLYcallablesitesFORANALYSIS.{sampling}.bedopssorted.bed"
	params:
		sort_bed = sort_bed_path
	shell:
		"{params.sort_bed} {input} > {output}"

rule intersect_callable_sites_mmul_analysis:
	input:
		expand(
			"callable_sites/{sample}.mmul.ONLYcallablesitesFORANALYSIS.{{sampling}}.bedopssorted.bed",
			sample=macaque_samples)
	output:
		"callable_sites/combined.macaque.mmul.INTERSECTIONcallablesites.{sampling}.bedopssorted.bed"
	params:
		bedops = bedops_path
	shell:
		"bedops -i {input} | bedops --merge - > {output}"

rule intersect_callable_sites_pcoq_analysis:
	input:
		expand(
			"callable_sites/{sample}.pcoq.ONLYcallablesitesFORANALYSIS.{{sampling}}.bedopssorted.bed",
			sample=sifaka_samples)
	output:
		"callable_sites/combined.sifaka.pcoq.INTERSECTIONcallablesites.{sampling}.bedopssorted.bed"
	params:
		bedops = bedops_path
	shell:
		"bedops -i {input} | bedops --merge - > {output}"

rule intersect_callable_sites_hg38_sifakas_analysis:
	input:
		expand(
			"callable_sites/{sample}.hg38.ONLYcallablesitesFORANALYSIS.{{sampling}}.bedopssorted.bed",
			sample=sifaka_samples)
	output:
		"callable_sites/combined.sifaka.hg38.INTERSECTIONcallablesites.{sampling}.bedopssorted.bed"
	params:
		bedops = bedops_path
	shell:
		"bedops -i {input} | bedops --merge - > {output}"

rule intersect_callable_sites_hg38_macaques_analysis:
	input:
		expand(
			"callable_sites/{sample}.hg38.ONLYcallablesitesFORANALYSIS.{{sampling}}.bedopssorted.bed",
			sample=macaque_samples)
	output:
		"callable_sites/combined.macaque.hg38.INTERSECTIONcallablesites.{sampling}.bedopssorted.bed"
	params:
		bedops = bedops_path
	shell:
		"bedops -i {input} | bedops --merge - > {output}"

rule sort_chrom_ascend_intersection:
	input:
		"callable_sites/combined.{species}.{genome}.INTERSECTIONcallablesites.{sampling}.bedopssorted.bed"
	output:
		"callable_sites/combined.{species}.{genome}.INTERSECTIONcallablesites.{sampling}.CHROMsorted.bed"
	shell:
		"sort -k1,1V -k2,2n {input} > {output}"
