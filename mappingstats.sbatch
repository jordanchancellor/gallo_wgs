#!/bin/bash
#SBATCH --mem=50gb
#SBATCH --time=1-0
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=12
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH -J mappingstats
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

# usage: sbatch --array=1-3 scripts/mappingstats.sbatch

# define global variables
wd="/project/noujdine_61/jchancel/gallo_oa_popgen_pipeline/"
samples_file="${wd}/samples_file.txt"
sample_id="$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)"
outdir="${wd}/var-calling/mpileups/"
bam="${outdir}/${sample_id}.bam"

# activate conda module and environments

module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
module load gcc/11.3.0
module load samtools

conda activate samtools

# create mapping stats file using samtools

echo "Creating mapping rate file with samtools flagstat for $bam."

samtools flagstat $bam > ${outdir}/${sample_id}.flagstat.txt

conda deactivate

echo "Done. Mapping stats files can be found in ${outdir}/${sample_id}.flagstat.txt."
