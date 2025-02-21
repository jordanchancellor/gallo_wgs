#!/bin/bash
#SBATCH --mem=24g
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=12
#SBATCH --partition=epyc-64
#SBATCH --error=outfiles/gatk2/%x_%a.err
#SBATCH --output=outfiles/gatk2/%x_%a.out
#SBATCH -J validate_sams
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

index="/scratch2/jchancel/gallo_oa_popgen/var-calling/reference_data/GCA_025277285.1_MgalMED_genomic.fna"
trimmedreads="/scratch2/jchancel/gallo_oa_popgen/gallo_trimmed/"
outdir="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/"
samples_file="/scratch2/jchancel/gallo_oa_popgen/bwa_samples.txt"
sample_id="$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)"
reads1="${trimmedreads}/${sample_id}/${sample_id}_R1.fq.gz"
reads2="${trimmedreads}/${sample_id}/${sample_id}_R2.fq.gz"
samfile="var-calling/alignment/bwa/${sample_id}.sam"
bamfile="var-calling/alignment/bwa/${sample_id}.sorted.bam"
cleanbam="var-calling/alignment/gatk2/${sample_id}.cleaned.sorted.bam"

#activate conda environment and load required modules
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate gatk
module load gcc/11.3.0
module load picard/2.26.2
module load gatk

# Run GATK4's ValidateSamFile in summary mode
gatk ValidateSamFile --INPUT ${cleanbam} \
--REFERENCE_SEQUENCE ${index} \
--MODE SUMMARY \
--OUTPUT ${outdir}/${sample_id}.validate.summary

# Run GATK4's ValidateSamFile in verbose mode
#gatk ValidateSamFile --INPUT ${cleanbam} \
#--REFERENCE_SEQUENCE ${index} \
#--MODE VERBOSE \
#--IGNORE_WARNINGS true \
#--OUTPUT ${outdir}/${sample_id}.validate2.verbose
