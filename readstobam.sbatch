#!/bin/bash
#SBATCH --mem=10GB
#SBATCH --time=10:00:00
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=12
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH -J bambam
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

# a pipeline to map trimmed fastq reads to a reference genome and produce cleaned, sorted bam files
# a bash and conda-dependent pipeline

# usage: sbatch --array=1-9 scripts/readstobam.sbatch

# define global variables and checkpoints

wd="/project/noujdine_61/jchancel/gallo_oa_popgen_pipeline/"
outdir="${wd}/var-calling/alignment/"
index="/project/noujdine_61/jchancel/indices/mussel/mytilus/galloprovinciales/MGAL_10.fa"
trimmedreads="${wd}/gallo_trimmed/"
samples_file="${wd}/bwa_samples.txt"
sample_id="$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)"
reads1="${trimmedreads}/${sample_id}/${sample_id}_R1.fq.gz"
reads2="${trimmedreads}/${sample_id}/${sample_id}_R2.fq.gz"
# sai1="${outdir}/${sample_id}_R1.sai"
# sai2="{outdir}/${sample_id}_R2.sai"
header="$(zcat ${reads1} | head -n 1)"
id="$(echo ${header} | awk '{print $1}' | sed 's/@//g' | cut -f 1-4 -d ":" | cut -f 3,4 -d ":" | sed 's/:/./g')"
sm="${sample_id}"
lb="$(echo ${header} | awk '{print $2}' | cut -f 4 -d ":")"
samfile="${outdir}/${sample_id}.sam"
bamfile="${outdir}/${sample_id}.sorted.bam"
calmdbam="${outdir}/${sample_id}.calmd.bam"
cleanbam="${outdir}/alignment/${sample_id}.cleaned.sorted.bam"

# activate conda module and environments
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate bwa

# uncomment code chunks below to use desired bwa algorithm

# align reads using bwa aln  
# echo "Running bwa aln on read1 files associated with ${reads1}."

# cmd1="bwa aln $index $reads1 -t $SLURM_CPUS_PER_TASK"

# echo "cmd1"

# $cmd1 > "${outdir}/${sample_id}_R1.sai"

# echo "Running bwa aln on read2 files associated with ${reads2}."

# cmd2="bwa aln $index $reads2 -t $SLURM_CPUS_PER_TASK"

# echo "$cmd2"

# $cmd2 > "${outdir}/${sample_id}_R2.sai"

# # convert sai to sam file using bwa sampe for paired-end reads

# echo "Read Group @RG\tID:${id}\tSM:${sm}\tLB:${lb}\tPL:ILLUMINA"

# echo "Running bwa sampe on 2 read files associated with ${reads1} and ${reads2}."

# cmdsampe="bwa sampe -r $(echo "@RG\tID:$id\tSM:$sm\tLB:$lb\tPL:ILLUMINA") $index $sai1 $sai2 $reads1 $reads2"

# echo "$cmdsampe"

# $cmdsampe > "${outdir}/${sample_id}.sam"

# convert sam files to sorted bam files using samtools

# align reads using bwa mem
echo "Running bwa mem on files associated with ${reads1} and ${reads2}"

bwacmd="bwa mem -t $SLURM_CPUS_PER_TASK -R $(echo "@RG\tID:$id\tSM:$sm\tLB:$lb\tPL:ILLUMINA") $index $reads1 $reads2"

echo "$bwacmd" 

$bwacmd > "${outdir}/${sample_id}.sam"

conda deactivate
conda activate samtools

echo "Converting sam files to sorted bam files with samtools."

samfile="${outdir}/${sample_id}.sam"

cat ${samfile} | samtools view -bSh | samtools sort -o ${outdir}/${sample_id}.sorted.bam

echo "Creating samtools index file from ${bamfile}"

samtools index -b ${bamfile} 

echo "Running idxstats on sorted, indexed bam file from ${bamfile}"

samtools idxstats ${bamfile} | tee ${outdir}/${sample_id}.idxstats.txt

echo "Creating mapping rate file with samtools flagstat from ${bamfile}"

samtools flagstat ${bamfile} > ${outdir}/${sample_id}.flagstat.txt

echo "Done."

conda deactivate

# validate sam files 

conda activate gatk
module load gcc/11.3.0
module load picard/2.26.2
module load gatk

# Run GATK4's ValidateSamFile in summary mode
gatk ValidateSamFile --INPUT ${bamfile} \
--REFERENCE_SEQUENCE ${index} \
--MODE SUMMARY \
--OUTPUT ${outdir}/${sample_id}.validate.summary

conda deactivate

# clean bam files if necessary 

if grep -q "ERROR:INVALID_TAG_NM" ${outdir}/${sample_id}.validate.summary ; then

    conda activate samtools
    module load gcc/11.3.0
    module load samtools

    samtools calmd -bAr ${bamfile} ${index} > ${outdir}/${sample_id}.calmd.bam
else
    echo "ERROR:INVALID_TAG_NM not found."
fi

conda deactivate

if grep -qe "#ERROR:MISMATCH_FLAG_MATE_NEG_STRAND" -e "ERROR:MISMATCH_FLAG_MATE_UNMAPPED" ${outdir}/${sample_id}.validate.summary ; then
    conda activate gatk
    module load picard/2.26.2
    module load gatk

    picard CleanSam --INPUT ${calmdbam} \
    --OUTPUT ${outdir}/${sample_id}.cleaned.sorted.bam
else
    echo "#ERROR:MISMATCH_FLAG_MATE_NEG_STRAND and ERROR:MISMATCH_FLAG_MATE_UNMAPPED not found."
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
    --OUTPUT ${outdir}/${sample_id}.cleanbam.revalidate.summary
else
    echo "${cleanbam} does not exist."
    if [ -f "${calmdbam}" ]; then
        echo "${calmdbam} exists. Re-validating bam file"
        # re-evaluate files 

        conda activate gatk
        module load gcc/11.3.0
        module load picard/2.26.2
        module load gatk

         # Run GATK4's ValidateSamFile in summary mode
        gatk ValidateSamFile --INPUT ${calmdbam} \
        --REFERENCE_SEQUENCE ${index} \
        --MODE SUMMARY \
        --OUTPUT ${outdir}/${sample_id}.calmdbam.revalidate.summary

    else
        echo "${calmdbam} does not exist."
    fi
fi

conda deactivate

echo "Pipeline complete. Trimmed fastq reads sucessfully mapped to reference genome and bam files cleaned and sorted."
echo "Disco!"