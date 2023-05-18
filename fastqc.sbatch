#!/bin/bash
#SBATCH --mem=5gb
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=12

# Load conda enviornments and fastqc
module load conda
conda activate fastqc
fastqc --version
if [[ $? -ne 0 ]]
then
	printf "Error - check FASTQC installation.\n"
	exit 1
fi

for i in $samples_file
do
fastqc $i
done

printf "Running FASTQC on read files associated \
with $i"
