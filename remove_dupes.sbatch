#!/bin/bash
#SBATCH --mem=50g
#SBATCH --time=42:00:00
#SBATCH --cpus-per-task=12
#SBATCH --partition=gpu
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH -J gatk_remove_dupes
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

# a pipeline to remove duplicate reads from sorted and cleaned bam files with gatk, call snps with haplotype caller, and filter variants
# a bash and conda-dependent pipeline

# usage: sbatch --array=1-9 scripts/remove_dupes.sbatch 

# define global variables and checkpoints

wd="/project/noujdine_61/jchancel/gallo_oa_popgen_pipeline/"
outdir="${wd}/var-calling/alignment/"
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
# samfile="${outdir}/${sample_id}.sam"
# bamfile="${outdir}/${sample_id}.sorted.bam"
# calmdbam="${outdir}/${sample_id}.calmd.bam"
cleanbam="${outdir}/${sample_id}.calmd.bam"
dedupbamfile="${outdir}/${sample_id}_marked_duplicates.bam"
# targetint="${outdir}/${sample_id}.realignertargetcreator.intervals"
rawvcf="${outdir}/gallo_OA.vcf.gz"
rawsnps="${outdir}/raw_snps_recal.vcf.gz"
snps="${outdir}/gallo_snps_passed_sites.vcf.gz"


#activate conda environment and load required modules
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh

# # build reference dictionary

# conda activate samtools

# echo "Building reference dictionary, fasta index, and bam index."

# samtools faidx ${index}

# conda deactivate

# create sequence dictionary 

conda activate gatk
module load gcc/11.3.0
module load picard/2.26.2
module load gatk

picard CreateSequenceDictionary -R ${index}

echo "Marking and removing duplicate reads from sorted bam file ${cleanbam} ith Picard."

picard MarkDuplicates \
-I ${cleanbam} \
-O ${outdir}/${sample_id}_marked_duplicates.bam \
-M ${outdir}/${sample_id}_marked_dup_metrics.txt \
--CREATE_INDEX true \
--REMOVE_DUPLICATES true \
--CREATE_MD5_FILE true

# run indel realignment

# echo "Running indel realignment on ${dedupbamfile}."

# java -jar /project/noujdine_61/rachelh/env/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
# -T RealignerTargetCreator \
# -R  ${index} \
# -I ${dedupbamfile} \
# -o ${outdir}/${sample_id}.realignertargetcreator.intervals

# java -jar /project/noujdine_61/rachelh/env/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
# -T IndelRealigner \
# -R ${index} \
# -targetIntervals ${targetint} \
# -I ${dedupbamfile} \
# -o ${outdir}/${sample_id}.indelrealigner.bam

# fix mate information

echo "Fixing mate information on ${dedupbamfile}."

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

gatk --java-options HaplotypeCaller \
-R ${index} \
-I ${outdir}/${sample_id}.fixed_mate.bam \
-O ${outdir}/${sample_id}.vcf.gz \
-ERC GVCF

# combine gvcfs from haplotype caller into single gvcf file with gatk

echo "Combining gvcfs from haplotype caller into single gvcf file with gatk."

gatk --java-options "-Xmx60g" CombineGVCFs \
-R ${index} \
-V ${outdir}/MG22_D2_1.g.vcf.gz \
-V ${outdir}/MG22_D2_2.g.vcf.gz \
-V ${outdir}/MG22_D2_3.g.vcf.gz \
-V ${outdir}/MG22_D23_C_1.g.vcf.gz \
-V ${outdir}/MG22_D23_C_2.g.vcf.gz \
-V ${outdir}/MG22_D23_C_3.g.vcf.gz \
-V ${outdir}/MG22_D25_T_1.g.vcf.gz \
-V ${outdir}/MG22_D25_T_2.g.vcf.gz \
-V ${outdir}/MG22_D25_T_3.g.vcf.gz \
-O ${outdir}/gallo_cohort.g.vcf.gz

# genotype gvcf file with gatk to obtain single, cohort vcf file

echo "Genotyping ${outdir}/gallo_cohort.g.vcf.gz with gatk to obtain single, cohort vcf file."

gatk --java-options "-Xmx60g" GenotypeGVCFs \
-R ${index} \
-V ${outdir}/gallo_cohort.g.vcf.gz \
-O ${outdir}/gallo_OA.vcf.gz

# extract snps with gatk

echo "Extracting SNPs from ${rawvcf} with gatk SelectVariants."

gatk SelectVariants \
        -R ${index} \
        -V ${rawvcf} \
        -select-type-to-include SNP \
        -O ${outdir}/raw_snps_recal.vcf.gz

echo "Filtering SNPs"

gatk VariantFiltration \
-V ${rawsnps} \
--filter-expression "!vc.hasAttribute('DP')" \
--filter-name "noCoverage" \
--filter-expression "vc.hasAttribute('DP') && DP < MINDEPTH" \
--filter-name "MinCov" \
--filter-expression "vc.hasAttribute('DP') && DP > MAXDEPTH" \
--filter-name "MaxCov" \
--filter-expression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || (vc.hasAttribute('QD') && QD < 2.0) " \
--filter-name "badSeq" \
--filter-expression "(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))" \
--filter-name "badStrand" \
--filter-expression "vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))" \
--filter-name "badMap" \
-O ${outdir}/filtered_snps.vcf.gz

#Extract only filter passed sites into hard filtered vcf

echo "Hard filtering ${outdir}/filtered_snps.vcf.gz with gatk."

gatk SelectVariants \
-V ${outdir}/filtered_snps.vcf.gz \
--exclude-filtered true \
-O ${outdir}/gallo_snps_passed_sites.vcf.gz

echo "Running gatk's validate variants on ${snps}."

# validate variants with gatk

gatk ValidateVariants \
   --validation-type-to-exclude ALLELES \
   -R ${index} \
   -V ${snps} > \
   ${outdir}/snps.validate.summary 2>&1

echo "Outputting sample-specific vcf files."

#output sample-specific vcfs

gatk SelectVariants \
-V ${snps} \
-sn MG22_D23_C_1 \
-sn MG22_D23_C_2 \
-sn MG22_D23_C_3 \
-O ${outdir}/MG22_D23_C_snps_passedfilteredsites.vcf.gz

gatk SelectVariants \
-V ${snps} \
-sn MG22_D25_T_1 \
-sn MG22_D25_T_2 \
-sn MG22_D25_T_3 \
-O ${outdir}/MG22_D25_T_snps_passedfilteredsites.vcf.gz

gatk SelectVariants \
-V ${snps} \
-sn MG22_D2_1 \
-sn MG22_D2_2 \
-sn MG22_D2_3 \
-O ${outdir}/MG22_D2_snps_passedfilteredsites.vcf.gz

conda deactivate

echo "Disco!"
