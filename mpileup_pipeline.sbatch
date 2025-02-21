#!/bin/bash
#SBATCH --mem=300gb
#SBATCH --time=2-0
#SBATCH --partition=largemem
#SBATCH --nodes=1
#SBATCH --nodelist=a01-10,a02-10,a04-10
#SBATCH --cpus-per-task=12
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH -J getbamfiles
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

# a pipeline to map trimmed fastq reads to a reference genome, produce cleaned, sorted, quality filtered bam files, and call snps using both gatk HaplotypeCaller and samtools mpileup
# a bash and conda-dependent pipeline

# usage: sbatch --array=1-3 scripts/mpileup_pipeline.sbatch

# define global variables
wd="/project/noujdine_61/jchancel/gallo_oa_popgen_pipeline/"
outdir="${wd}/var-calling/mpileups"
index="/project/noujdine_61/jchancel/indices/mussel/mytilus/galloprovinciales/JAWDJN01.fa"
trimmedreads="${wd}/gallo_trimmed/"
samples_file="${wd}/samples_file.txt"
sample_id="$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)"
reads1="${trimmedreads}/${sample_id}/${sample_id}_R1.fq"
reads2="${trimmedreads}/${sample_id}/${sample_id}_R2.fq"
header="$(zcat ${reads1} | head -n 1)"
id="$(echo ${header} | awk '{print $1}' | sed 's/@//g' | cut -f 1-4 -d ":" | cut -f 3,4 -d ":" | sed 's/:/./g')"
sm="${sample_id}"
lb="$(echo ${header} | awk '{print $2}' | cut -f 4 -d ":")"
# mpileupworkflow variables
samfile="${outdir}/${sample_id}.sam"
bamfile="${outdir}/${sample_id}.bam"
sortedbamfile="${outdir}/${sample_id}.sorted.bam"
dupesremovedbamfile="${outdir}/${sample_id}_duplicates_removed.sorted.bam"
dupesremovedqfbamfile="${outdir}/${sample_id}_duplicates_removed.qf.sorted.bam"

# activate conda module and environments

module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
module load gcc/11.3.0
module load picard/2.26.2
module load gatk
module load samtools

# conda activate bwa

# # align reads using bwa mem
# echo "Running bwa mem on files associated with ${reads1} and ${reads2}"

# bwacmd="bwa mem -t $SLURM_CPUS_PER_TASK -R $(echo "@RG\tID:$id\tSM:$sm\tLB:$lb\tPL:ILLUMINA") $index $reads1 $reads2"

# echo "$bwacmd" 

# $bwacmd > "${outdir}/${sample_id}.sam"

# conda deactivate

# # prepare bam file for snp calling

# conda activate samtools

# # convert sam file from bwa alignment to bam file using samtools

# echo "Creating bam file from ${samfile}"

# samtools view -Sb ${samfile} > ${outdir}/${sample_id}.bam

# conda deactivate
conda activate gatk

# sort bam file with picard

echo "Sorting ${bamfile} with Picard"

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
conda activate samtools

# remove low aquality alignments from ambiguous mapping with samtools

echo "Removing low aquality alignments from ${dupesremovedbamfile} with Samtools"

samtools view -b -f 0x0002 -F 0x0004 -F 0x0008 -q 20 ${dupesremovedbamfile} > ${outdir}/${sample_id}_duplicates_removed.qf.sorted.bam

# create mapping stats file

echo "Creating mapping rate file with samtools flagstat from ${dupesremovedqfbamfile}"

samtools flagstat ${dupesremovedqfbamfile} > ${outdir}/${sample_id}.picard.flagstat.txt

conda deactivate
conda activate gatk

# validate bam files

echo "Validating ${dupesremovedqfbamfile} with GATK4."

# Run GATK4's ValidateSamFile in summary mode
gatk ValidateSamFile --INPUT ${dupesremovedqfbamfile} \
--REFERENCE_SEQUENCE ${index} \
--MODE SUMMARY \
--OUTPUT ${outdir}/${sample_id}_duplicates_removed.qf.sorted.validate.summary

conda deactivate

echo "DISCO!"
