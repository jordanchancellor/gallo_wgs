#!/bin/bash
#SBATCH --mem=50GB
#SBATCH --time=10:00:00
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=12
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH -J fixmateinfo
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

# fix mate information on bam files with picard
# bash conda-dependent

# usage: sbatch --array=1-9 scripts/cleanbamfiles.sbatch

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
# samfile="${wd}/var-calling/alignment/${sample_id}.sam"
# bamfile="${outdir}/${sample_id}.bam"
# sortedbamfile="${outdir}/${sample_id}.sorted.bam"
# dupesremovedbamfile="${outdir}/${sample_id}_duplicates_removed.sorted.bam"
dupesremovedqfbamfile="${outdir}/${sample_id}_duplicates_removed.qf.sorted.bam"
fixedmatesbam="${outdir}/${sample_id}_duplicates_removed.qf.matesfixed.sorted.bam"
cleanbam="${outdir}/${sample_id}_duplicates_removed.qf.matesfixed.sorted.calmd.bam"

# activate conda module and environments
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate gatk
module load gcc/11.3.0
module load picard/2.26.2
module load gatk

# fix mate information

picard FixMateInformation \
       --INPUT ${dupesremovedqfbamfile} \
       --OUTPUT ${outdir}/${sample_id}_duplicates_removed.qf.matesfixed.sorted.bam \
       --VALIDATION_STRINGENCY LENIENT


# run gatk ValidateSamFile

gatk ValidateSamFile --INPUT ${fixedmatesbam} \
        --REFERENCE_SEQUENCE ${index} \
        --MODE SUMMARY \
        --OUTPUT ${outdir}/${sample_id}_duplicates_removed.qf.matesfixed.sorted.validated.summary

# clean bam files if necessary 

if grep -q "ERROR:INVALID_TAG_NM" ${outdir}/${sample_id}_duplicates_removed.qf.matesfixed.sorted.validated.summary ; then

    conda activate samtools
    module load gcc/11.3.0
    module load samtools

    samtools calmd -bAr ${fixedmatesbam} ${index} > ${outdir}/${sample_id}_duplicates_removed.qf.matesfixed.sorted.calmd.bam
else
    echo "ERROR:INVALID_TAG_NM not found."
fi

conda deactivate

if [ -f "${cleanbam}" ]; then
    echo "${cleanbam} exists. Re-validating bam file"

    # re-evaluate files 

    conda activate gatk
    module load gcc/11.3.0
    module load picard/2.26.2
    module load gatk

    # Run GATK4's ValidateSamFile in summary mode
    gatk ValidateSamFile --INPUT ${cleanbam} \
    --REFERENCE_SEQUENCE ${index} \
    --MODE SUMMARY \
    --OUTPUT ${outdir}/${sample_id}_duplicates_removed.qf.matesfixed.sorted.calmd.validate.summary
else
    echo "${cleanbam} does not exist."
fi

conda deactivate
