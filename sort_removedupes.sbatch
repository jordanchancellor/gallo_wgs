#!/bin/bash
#SBATCH --mem=36GB
#SBATCH --time=05:00:00
#SBATCH --partition=epyc-64
#SBATCH --cpus-per-task=16
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH -J samtools_sort
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

# convert sam to bam file, sort, mark and remove duplicates, remove low quality alignments
# bash conda-dependent

# usage: sbatch --array=1-9 scripts/sort_removedupes.sbatch

# define global variables and checkpoints

wd="/project/noujdine_61/jchancel/gallo_oa_popgen_pipeline/"
outdir="${wd}/var-calling/alignment/picard/"
index="/project/noujdine_61/jchancel/indices/mussel/mytilus/galloprovinciales/MGAL_10.fa"
# trimmedreads="${wd}/gallo_trimmed/"
samples_file="${wd}/bwa_samples.txt"
sample_id="$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)"
# reads1="${trimmedreads}/${sample_id}/${sample_id}_R1.fq.gz"
# reads2="${trimmedreads}/${sample_id}/${sample_id}_R2.fq.gz"
# sai1="${outdir}/${sample_id}_R1.sai"
# sai2="{outdir}/${sample_id}_R2.sai"
# header="$(zcat ${reads1} | head -n 1)"
# id="$(echo ${header} | awk '{print $1}' | sed 's/@//g' | cut -f 1-4 -d ":" | cut -f 3,4 -d ":" | sed 's/:/./g')"
# sm="${sample_id}"
# lb="$(echo ${header} | awk '{print $2}' | cut -f 4 -d ":")"
samfile="${wd}/var-calling/alignment/${sample_id}.sam"
bamfile="${outdir}/${sample_id}.bam"
sortedbamfile="${outdir}/${sample_id}.sorted.bam"
dupesremovedbamfile="${outdir}/${sample_id}_duplicates_removed.sorted.bam"
dupesremovedqfbamfile="${outdir}/${sample_id}_duplicates_removed.qf.sorted.bam"
# calmdbam="${outdir}/${sample_id}.calmd.bam"
# cleanbam="${outdir}/alignment/${sample_id}.cleaned.sorted.bam"

# activate conda module and environments
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate samtools

# convert sam file from bwa alignment to bam file using samtools

echo "Creating bam file from ${samfile}"

samtools view -Sb ${samfile} > ${outdir}/${sample_id}.bam

conda deactivate

# sort bam file with picard

echo "Sorting ${bamfile} with Picard"

conda activate gatk
module load gcc/11.3.0
module load picard/2.26.2

picard SortSam -I ${bamfile} \
-O ${outdir}/${sample_id}.sorted.bam \
--VALIDATION_STRINGENCY SILENT \
--SORT_ORDER coordinate

# remove duplicates from sorted bam file with picard 

echo "Marking and removing duplicates from ${sortedbamfile} with Picard"

picard MarkDuplicates \
-I ${sortedbamfile} \
-O ${outdir}/${sample_id}_duplicates_removed.sorted.bam \
-M ${outdir}/${sample_id}_marked_duplicate_bam_metrics.txt \
--VALIDATION_STRINGENCY SILENT \
--CREATE_INDEX true \
--REMOVE_DUPLICATES true \
--CREATE_MD5_FILE true

conda deactivate

# remove low aquality alignments from ambiguous mapping with samtools

conda activate samtools

echo "Removing low aquality alignments from ${dupesremovedbamfile} with Picard"

samtools view -b -f 0x0002 -F 0x0004 -F 0x0008 -q 20 ${dupesremovedbamfile} > ${outdir}/${sample_id}_duplicates_removed.qf.sorted.bam

# create mapping stats file

echo "Creating mapping rate file with samtools flagstat from ${dupesremovedqfbamfile}"

samtools flagstat ${dupesremovedqfbamfile} > ${outdir}/${sample_id}.picard.flagstat.txt

conda deactivate

# validate bam files

echo "Validating ${dupesremovedqfbamfile} with GATK4."

conda activate gatk
module load gcc/11.3.0
module load picard/2.26.2
module load gatk

# Run GATK4's ValidateSamFile in summary mode
gatk ValidateSamFile --INPUT ${dupesremovedqfbamfile} \
--REFERENCE_SEQUENCE ${index} \
--MODE SUMMARY \
--OUTPUT ${outdir}/${sample_id}_duplicates_removed.qf.sorted.validate.summary

conda deactivate

echo "Disco!"
