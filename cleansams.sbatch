#!/bin/bash
#SBATCH --mem=24g
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=12
#SBATCH --partition=epyc-64
#SBATCH --error=outfiles/gatk2/%x_%a.err
#SBATCH --output=outfiles/gatk2/%x_%a.out
#SBATCH -J cleansams
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

index="/scratch2/jchancel/gallo_oa_popgen/var-calling/reference_data/GCA_025277285.1_MgalMED_genomic.fna"
outdir="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/"
samples_file="/scratch2/jchancel/gallo_oa_popgen/bwa_samples.txt"
sample_id="$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)"
samfile="var-calling/alignment/bwa/${sample_id}.sam"
bamfile="var-calling/alignment/bwa/${sample_id}.sorted.bam"
calmdbam="var-calling/alignment/gatk2/${sample_id}.calmd.bam"

#load conda environments and modules
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate samtools
module load gcc/11.3.0
module load samtools

#clean bam files
samtools calmd -bAr ${bamfile} ${index} > ${outdir}/${sample_id}.calmd.bam

conda deactivate
conda activate gatk
module load picard/2.26.2
module load gatk

picard CleanSam --INPUT ${calmdbam} \
--OUTPUT ${outdir}/${sample_id}.cleaned.sorted.bam

