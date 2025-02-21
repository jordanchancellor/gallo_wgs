#!/bin/bash
#SBATCH --mem=32gb
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=12
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out
#SBATCH -J bwa_index

refdir="/project/noujdine_61/jchancel/indices/mussel/mytilus/galloprovinciales/"
index="/project/noujdine_61/jchancel/indices/mussel/mytilus/galloprovinciales/JAWDJN01.fa"

#activate conda environments
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate bwa

printf "Running bwa index.\n"

#create bwa index
cd $refdir || exit
bwa index -a bwtsw $index
