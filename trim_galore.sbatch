#!/bin/bash
#SBATCH --mem=5gb
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=12
#SBATCH -J trimgalore
#SBATCH -o %x_%a.out

samples_dir='01.RawData/*'
trimmed_dir='gallo_trimmed'
samples_file='samples_file.txt'
sample_id=`cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p`
reads1="${samples_dir}/${sample_id}_1.fq.gz"
reads2="${samples_dir}/${sample_id}_2.fq.gz"

# Load trim_galore
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate trimgalore
trim_galore --version
if [[ $? -ne 0 ]]
then
 	printf "Error - check Trim Galore! installation.\n"
fi

printf "Running Trim Galore! on 2 read files associated \
with ${sample_id}: \n\
$reads1 and ${reads2}.\n\n"

trim_galore --paired --phred33 --output_dir $trimmed_dir \
--length 36 --stringency 1 -e 0.1 \
--fastqc_args "-t ${SLURM_CPUS_PER_TASK}" $reads1 $reads2

# Checkpoint
if [[ $? -eq 0 ]]
then
	if [[ ! -d $trimmed_dir ]]
	then
		mkdir $trimmed_dir
	fi

#	if [[ `ls ${outdir}/${sample_id}*val* | wc -l` == 6 ]] && \
#[[ -d $trimmed_dir ]]
#	then
#		printf "Trimming and FastQC run completed \
#for ${sample_id}.\n"
#		mv ${outdir}/${sample_id}*fq.gz $trimmed_dir
#		if [[ ! -d checkpoints ]]
#		then
#			mkdir checkpoints
#		fi
#		touch \
#checkpoints/${prefix}_${SLURM_ARRAY_TASK_ID}.checkpoint
#		exit 0
#	else
#		printf "Error - not all reads trimmed \
#for ${sample_id}.\n"
#		exit 1
#	fi
#else
#	printf "Error - non-zero exit code returned by Trim Galore!\n"
#	exit 1
#fi

