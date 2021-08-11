configfile: "config.yaml"

rule all:
	input:
		"plots/quals.svg"

include: "bwa_snakefile"


rule samtools sort:
	input:
		"mapped_reads/{sample}.bam"
	output:
		protected("sorted_reads/{sample}.bam")
	log:
		"logs/samtools/{sample}_sort.log"
	shell:
		"samtools sort -T sorted_reads/{wildcards.sample} "
		"-O bam {input} > {output} 2> {log}"

rule samtools_index:
	input:
		"sorted_reads/{sample}.bam"
	output:
		"sorted_reads/{sample}.bam.bai"
	log:
		"logs/samtools/{sample}_index.log"
	shell:
		"samtools index {input} 2> {log}"

rule bcftools_call:
	input:
		fa="data/genome.fa",
		bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
		bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
	output:
		"calls/all.vcf"
	params:
		p=config["mutationRate"]
	log:
		"logs/bcftools/all_samples.log"
	shell:
		"bcftools mpileup  -f {input.fa} {input.bam} | "
		"bcftools call -mv -P '{params.p}' - > {output} 2> {log}"

rule plot_quals:
	input:
		"calls/all.vcf"
	output:
		"plots/quals.svg"
	script:
		"scripts/plot-quals.py"
