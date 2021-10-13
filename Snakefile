# Author: Brian Wray
# Date: 13-Oct-2021
# ATACseq pipeline
# See README for details

configfile: "config.yaml"
curr_samples = config["samples"]

rule all:
	input:
		expand("analysis/mapped_reads/{sample}.bam",  sample = curr_samples),
		expand("data/fastqc/{sample}.fastqc.zip", sample = curr_samples)

rule fastqc:
	input:
		expand("data/reads/{sample}_R1.fq.gz", sample=curr_samples),
		expand("data/reads/{sample}_R2.fq.gz", sample=curr_samples)
	output:
		"data/fastqc/{sample}.fastqc.zip"
	threads:
		config["threads"]["fastqc"]
	params:
		out_dir = config["fastqc"]["outdir"]
	shell:
		"fastqc -o {params.out_dir} -t {threads} {input}"
	

rule trimGalore:
	input:
		expand("data/reads/{sample}_R1.fq.gz", sample=curr_samples),
		expand("data/reads/{sample}_R2.fq.gz", sample = curr_samples)
	output:
		"data/trimmed_reads/{sample}_R1_trimmed.fq.gz",
		"data/trimmed_reads/{sample}_R2_trimmed.fq.gz"
	params:
		stringency = config["trim_galore"]["stringency"],
		qscore = config["trim_galore"]["qScore"],
		trimErrorRate = config["trim_galore"]["trimErrorRate"],
		paired = config["trim_galore"]["paired"],
		length = config["trim_galore"]["trimMinLength"],
		out_dir = config["trim_galore"]["out_dir"]
	log:
		"logs/trim_galore/{sample}.log"
	threads:
		config["threads"]["trimGalore"]
	shell:
		'''
		trim_galore -suppress_warn -o {params.out_dir} \
		--stringency {params.stringency} \
		-e {params.trimErrorRate} \
		--length {params.length} \
		-q {params.qscore} \
		{params.paired} \
		{input} 2> {log}
		'''
		
rule bwa:
	input:
		"data/genome.fa",
		expand("data/trimmed_reads/{sample}_R1_trimmed.fq.gz", sample = curr_samples),
		expand("data/trimmed_reads/{sample}_R2_trimmed.fq.gz", sample = curr_samples)
	output:
		"analysis/mapped_reads/{sample}.bam"
	params:
                rg=r"@RG\tID:{sample}\tSM:{sample}"
	log:
                "logs/bwa_mem/{sample}.log"
	threads: 
		config["threads"]["bwa"] 
	shell:
                '''
		bwa mem {input} -R {params.rg} -t {threads} {input} | \
                samtools view -Sb - > {output} 2> {log}
		'''

