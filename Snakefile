# Author: Brian Wray
# Date: 13-Oct-2021
# ATACseq pipeline
# See README for details

configfile: "config.yaml"
curr_samples = config["samples"]

rule all:
	input:
		expand("analysis/mapped_reads/{sample}_filtered.bam.bai",  sample = curr_samples),
		expand("data/fastq_screen/{sample}", sample = curr_samples),
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
	
rule fastq_screen:
	input:
		expand("data/reads/{sample}_R1.fq.gz", sample = curr_samples)
	output:
		"data/fastq_screen/{sample}"
	params:
		aligner = config["fastq_screen"]["aligner"],
		conf = config["fastq_screen"]["conf"],
		ex_path = config["fastq_screen"]["ex_path"],
		outdir = config["fastq_screen"]["outdir"]
	threads:
		config["threads"]["fastq_screen"]
	shell:
		'''
		perl {params.ex_path} \
		--threads {threads} \
		--conf {params.conf} \
		--outdir {params.outdir} \
		--aligner {params.aligner} \
		{input}
		'''

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
		
rule bwa_sorted:
	input:
		"data/genome.fa",
		expand("data/trimmed_reads/{sample}_R1_trimmed.fq.gz", sample = curr_samples),
		expand("data/trimmed_reads/{sample}_R2_trimmed.fq.gz", sample = curr_samples)
	output:
		"analysis/mapped_reads/{sample}_sorted.bam"
	params:
                rg=r"@RG\tID:{sample}\tSM:{sample}"
	log:
                "logs/bwa_mem/{sample}.log"
	threads: 
		config["threads"]["bwa"] 
	shell:
                '''
		bwa mem {input} -R {params.rg} -t {threads} {input} | \
                samtools sort -Sb - > {output} 2> {log}
		'''

rule remove_pcr_duplicates:
	input:
		expand("analysis/mapped_reads/{sample}_sorted.bam", sample = curr_samples)
	output:
		bam="analysis/mapped_reads/{sample}_noDups.bam",
		dupes="analysis/duplicates/{sample}_dups.txt"
	log:
		"logs/pcr_dups/{sample}.log"
	params:
		picard_path = config["picard"]["path"]
	shell:
		'''
		{params.picard_path}/picard.jar MarkDuplicates \
		I = {input} \
		O = {output.bam} \
		M = {output.dupes} \
		REMOVE_DUPLICATES = true
		'''

rule index_no_dups:
	input:
		expand("analysis/mapped_reads/{sample}_noDups.bam", sample = curr_samples)
	output:
		"analysis/mapped_reads/{sample}_noDups.bam.bai"
	log:
		"logs/pcr_dups/{sample}_index.log"
	shell:
		'''
		samtools index {input}
		'''

rule remove_mito:
	input:
		expand("analysis/mapped_reads/{sample}_noDups.bam", sample = curr_samples),
		expand("analysis/mapped_reads/{sample}_noDups.bam.bai", sample = curr_samples)
	output:
		"analysis/mapped_reads/{sample}_noMito.bam"
	log:
		"logs/remove_mito/{sample}.log"
	shell:
		'''
		samtools idxstats {input} | \
		cut -f 1 | \
		grep -v chrM | \
		xargs samtools view -b {input} > {output} 2> {log}
		'''

rule index_remove_mito:
	input:
		expand("analysis/mapped_reads/{sample}_noMito.bam", sample = curr_samples)
	output:
		"analysis/mapped_reads/{sample}_noMito.bam.bai"
	log:
		"logs/remove_mito/{sample}_index.log"
	shell:
		'''
		samtools index {input}
		'''

rule get_proper_bam:
	input:
		expand("analysis/mapped_reads/{sample}_noMito.bam", sample = curr_samples),
		expand("analysis/mapped_reads/{sample}_noMito.bam.bai", sample = curr_samples),
	output:
		"analysis/mapped_reads/{sample}_filtered.bam"
	log:
		"logs/filtered_bam/{sample}.log"	
	shell:
		'''
		samtools view -h -q 30 {input} | \
		samtools view -h -b -F 1804 -f 2  > {output}
		'''

rule get_filtered_bam_idx:
	input:
		expand("analysis/mapped_reads/{sample}_filtered.bam", sample = curr_samples)
	output:
		"analysis/mapped_reads/{sample}_filtered.bam.bai"
	log:
		"logs/filtered_bam/{sample}_index.log"
	shell:
		'''
		samtools index {input}
		''' 
	

