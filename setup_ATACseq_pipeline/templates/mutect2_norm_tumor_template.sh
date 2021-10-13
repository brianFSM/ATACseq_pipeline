#!/bin/zsh
#SBATCH -A b1042             ## account (unchanged)
#SBATCH -p genomics          ## "-p" instead of "-q"
#SBATCH -J mutect2_hg38
#SBATCH --mail-type=FAIL,TIME_LIMIT_90
#SBATCH --mail-user=brian.wray@northwestern.edu
#SBATCH -o "%x.o%j"
#SBATCH -N 1                 ## number of nodes
#SBATCH -n 1                 ## number of cores
#SBATCH -t 30:00:00          ## walltime
#SBATCH --mem=100G
source /etc/zshrc

module purge all
module load gatk/4.1.0
module load java/jdk1.8.0_25

sample={{sample}}
reference_dir="/projects/b1012/xvault/REFERENCES/builds/{{reference}}"
bam_dir=analysis/bam_files
outfile=analysis/bam/${sample}_{{reference}}_somatic.vcf.gz

printf "Starting mutect2 on %s at " $sample
date

set -x
gatk Mutect2 \
	-R ${reference_dir}/genome.fa \
     {{tumor_files}} {{normal_files}} {{normal_samples}}	--germline-resource ${reference_dir}/bundle/somatic-hg38-af-only-gnomad.hg38.vcf.gz \
	-O $outfile 2> logs/${sample}_mutect2.log

set +x

printf "Finished running mutect2 on %s at " $sample
date
