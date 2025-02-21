#!/bin/bash
#SBATCH --mem=500GB
#SBATCH --time=7-0
#SBATCH --partition=largemem
#SBATCH --nodes=1
#SBATCH --nodelist=a02-10,a03-10,a04-10
#SBATCH --cpus-per-task=12
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH -J snpcalling
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

# a pipeline to map trimmed fastq reads to a reference genome, produce cleaned, sorted, quality filtered bam files, and call snps using both gatk HaplotypeCaller and samtools mpileup
# a bash and conda-dependent pipeline

# usage: sbatch --array=1-9 scripts/pipeline-mgallo-snp-calling.sbatch

# define global variables
wd="/project/noujdine_61/jchancel/gallo_oa_popgen_pipeline/"
outdir="${wd}/var-calling/alignment/gatk"
mpileupoutdir="${wd}/var-calling/alignment/mpileup"
popooldir="${wd}/var-calling/popoolation/"
popoolsoftware="/project/noujdine_61/jchancel/software/popoolation2_1201/"
index="/project/noujdine_61/jchancel/indices/mussel/mytilus/galloprovinciales/MGAL_10.fa"
trimmedreads="${wd}/gallo_trimmed/"
samples_file="${wd}/bwa_samples.txt"
sample_id="$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)"
reads1="${trimmedreads}/${sample_id}/${sample_id}_R1.fq.gz"
reads2="${trimmedreads}/${sample_id}/${sample_id}_R2.fq.gz"
header="$(zcat ${reads1} | head -n 1)"
id="$(echo ${header} | awk '{print $1}' | sed 's/@//g' | cut -f 1-4 -d ":" | cut -f 3,4 -d ":" | sed 's/:/./g')"
sm="${sample_id}"
lb="$(echo ${header} | awk '{print $2}' | cut -f 4 -d ":")"
# gatk workflow variables 
samfile="${outdir}/${sample_id}.sam"
sortedbamfile="${outdir}/${sample_id}.sorted.bam"
calmdbam="${outdir}/${sample_id}.calmd.sorted.bam"
cleanbam="${outdir}/${sample_id}.cleaned.sorted.bam"
dedupbamfile="${outdir}/${sample_id}_marked_duplicates.bam"
# mpileupworkflow variables
bamfile="${mpileupoutdir}/${sample_id}.bam"
mpileupsortedbamfile="${mpileupoutdir}/${sample_id}.sorted.bam"
dupesremovedbamfile="${mpileupoutdir}/${sample_id}_duplicates_removed.sorted.bam"
dupesremovedqfbamfile="${mpileupoutdir}/${sample_id}_duplicates_removed.qf.sorted.bam"
mpileupfile="${mpileupoutdir}/${sample_id}.mpileup"
mpileupsync="${popooldir}/${sample_id}.mpileup.sync"

# activate conda module and environments

module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
module load gcc/11.3.0
module load picard/2.26.2
module load gatk
module load samtools

conda activate bwa

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

echo "Creating samtools index file from ${sortedbamfile}"

samtools index -b ${sortedbamfile} 

echo "Running idxstats on sorted, indexed bam file from ${sortedbamfile}"

samtools idxstats ${sortedbamfile} | tee ${outdir}/${sample_id}.idxstats.txt

echo "Creating mapping rate file with samtools flagstat from ${bamfile}"

samtools flagstat ${sortedbamfile} > ${outdir}/${sample_id}.flagstat.txt

conda deactivate
conda activate gatk

# Run GATK4's ValidateSamFile in summary mode
gatk ValidateSamFile --INPUT ${sortedbamfile} \
--REFERENCE_SEQUENCE ${index} \
--MODE SUMMARY \
--OUTPUT ${outdir}/${sample_id}.validate.summary

conda deactivate

# call snps using GATK4 HaplotypeCaller workflow

# clean bam files if necessary 

if grep -q "ERROR:INVALID_TAG_NM" ${outdir}/${sample_id}.validate.summary ; then

    conda activate samtools

    samtools calmd -bAr ${sortedbamfile} ${index} > ${outdir}/${sample_id}.calmd.sorted.bam
else
    echo "ERROR:INVALID_TAG_NM not found."
fi

conda deactivate

if grep -qe "#ERROR:MISMATCH_FLAG_MATE_NEG_STRAND" -e "ERROR:MISMATCH_FLAG_MATE_UNMAPPED" ${outdir}/${sample_id}.validate.summary ; then
    
    conda activate gatk

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

         # Run GATK4's ValidateSamFile in summary mode
        gatk ValidateSamFile --INPUT ${calmdbam} \
        --REFERENCE_SEQUENCE ${index} \
        --MODE SUMMARY \
        --OUTPUT ${outdir}/${sample_id}.calmdbam.revalidate.summary

    else
        echo "${calmdbam} does not exist."
    fi
fi

picard CreateSequenceDictionary -R ${index}

if [ -f "${cleanbam}" ]; then
    echo "Marking and removing duplicate reads from sorted bam file ${cleanbam} with Picard."

    picard MarkDuplicates \
    -I ${cleanbam} \
    -O ${outdir}/${sample_id}_marked_duplicates.bam \
    -M ${outdir}/${sample_id}_marked_dup_metrics.txt \
    --CREATE_INDEX true \
    --REMOVE_DUPLICATES true \
    --CREATE_MD5_FILE true
else
    echo "${cleanbam} does not exist."
    if [ -f "${calmdbam}" ]; then
    echo "${calmdbam} exists. Marking and removing duplicate reads from sorted bam file ${calmdbam} ith Picard."

    picard MarkDuplicates \
    -I ${calmdbam} \
    -O ${outdir}/${sample_id}_marked_duplicates.bam \
    -M ${outdir}/${sample_id}_marked_dup_metrics.txt \
    --CREATE_INDEX true \
    --REMOVE_DUPLICATES true \
    --CREATE_MD5_FILE true
    fi
fi

picard FixMateInformation \
-I ${dedupbamfile} \
-O ${outdir}/${sample_id}.fixed_mate.bam \
--ADD_MATE_CIGAR true

echo "Indexing ${outdir}/${sample_id}.fixed_mate.bam using samtools."

conda deactivate 
conda activate samtools

samtools index ${outdir}/${sample_id}.fixed_mate.bam

conda deactivate
conda activate gatk

# call snps with gatk haplotype caller

echo "Calling variants using HaplotypeCaller."

gatk --java-options "-Xmx480g" HaplotypeCaller \
-R ${index} \
-I ${outdir}/${sample_id}.fixed_mate.bam \
-O ${outdir}/${sample_id}.g.vcf.gz \
-ERC GVCF

# call snps using samtools mpileup; popoolation workflow

conda deactivate
conda activate samtools

# convert sam file from bwa alignment to bam file using samtools

echo "Creating bam file from ${samfile}"

samtools view -Sb ${samfile} > ${mpileupoutdir}/${sample_id}.bam

conda deactivate
conda activate gatk

# sort bam file with picard

echo "Sorting ${bamfile} with Picard"

picard SortSam -I ${bamfile} \
-O ${mpileupoutdir}/${sample_id}.sorted.bam \
--VALIDATION_STRINGENCY SILENT \
--SORT_ORDER coordinate

# remove duplicates from sorted bam file with picard 

echo "Marking and removing duplicates from ${mpileupsortedbamfile} with Picard"

picard MarkDuplicates \
-I ${mpileupsortedbamfile} \
-O ${mpileupoutdir}/${sample_id}_duplicates_removed.sorted.bam \
-M ${mpileupoutdir}/${sample_id}_marked_duplicate_bam_metrics.txt \
--VALIDATION_STRINGENCY SILENT \
--CREATE_INDEX true \
--REMOVE_DUPLICATES true \
--CREATE_MD5_FILE true

conda deactivate
conda activate samtools

# remove low aquality alignments from ambiguous mapping with samtools

echo "Removing low aquality alignments from ${dupesremovedbamfile} with Picard"

samtools view -b -f 0x0002 -F 0x0004 -F 0x0008 -q 20 ${dupesremovedbamfile} > ${mpileupoutdir}/${sample_id}_duplicates_removed.qf.sorted.bam

# create mapping stats file

echo "Creating mapping rate file with samtools flagstat from ${dupesremovedqfbamfile}"

samtools flagstat ${dupesremovedqfbamfile} > ${mpileupoutdir}/${sample_id}.picard.flagstat.txt

conda deactivate
conda activate gatk

# validate bam files

echo "Validating ${dupesremovedqfbamfile} with GATK4."

# Run GATK4's ValidateSamFile in summary mode
gatk ValidateSamFile --INPUT ${dupesremovedqfbamfile} \
--REFERENCE_SEQUENCE ${index} \
--MODE SUMMARY \
--OUTPUT ${mpileupoutdir}/${sample_id}_duplicates_removed.qf.sorted.validate.summary

conda deactivate
conda activate samtools

# compile snps by creating samtools mpileup file

echo "Creating mpileup file from ${dupesremovedqfbamfile} with samtools."

samtools mpileup -B -Q 0 -f ${index} ${dupesremovedqfbamfile} > ${mpileupoutdir}/${sample_id}.mpileup

conda deactivate

# synchronize file with popoolation 

echo "Syncing ${mpileupfile} with popoolation."

perl ${popoolsoftware}/mpileup2sync.pl --fastq-type illumina --min-qual 20 --input ${mpileupfile} --output ${popooldir}/${sample_id}.mpileup.sync

# subsample to uniform coverage with popoolation

# echo "Subsampling ${mpileupsync} to uniform coverage with popoolation."

# perl ${popoolsoftware}/subsample-synchronized.pl \
# --input ${mpileupsync} \
# --output ${popooldir}/${sample_id}_ss10_idf.mpileup.sync \
# --target-coverage 10 \
# --max-coverage '2%' \
# --method withoutreplace

echo "DISCO!"
