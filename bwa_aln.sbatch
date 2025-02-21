#!/bin/bash
#SBATCH --mem=50g
#SBATCH --time=1-0
#SBATCH --cpus-per-task=12
#SBATCH --error=outfiles/bwa/%x_%a.err
#SBATCH --output=outfiles/bwa/%x_%a.out
#SBATCH -J bwa_aln
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

index="/scratch2/jchancel/gallo_oa_popgen/var-calling/reference_data/GCA_025277285.1_MgalMED_genomic.fna"
trimmedreads="/scratch2/jchancel/gallo_oa_popgen/gallo_trimmed/"
outdir="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/bwa/"
samples_file="/scratch2/jchancel/gallo_oa_popgen/bwa_samples.txt"
sample_id="$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)"
reads1="${trimmedreads}/${sample_id}/${sample_id}_R1.fq.gz"
reads2="${trimmedreads}/${sample_id}/${sample_id}_R2.fq.gz"


#Load bwa conda environment
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate bwa

echo "Running bwa aln on read1 files associated with ${reads1}."

cmd1="bwa aln $index $reads1 -t $SLURM_CPUS_PER_TASK"

echo "cmd1"

$cmd1 > "var-calling/alignment/bwa/${sample_id}_R1.sai"

echo "Running bwa aln on read2 files associated with ${reads2}."

cmd2="bwa aln $index $reads2 -t $SLURM_CPUS_PER_TASK"

echo "$cmd2"

$cmd2 > "var-calling/alignment/bwa/${sample_id}_R2.sai"
