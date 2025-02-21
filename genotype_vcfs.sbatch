#!/bin/bash
#SBATCH --mem=180GB
#SBATCH --time=1-0
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --error=%x.err
#SBATCH --output=%x.out
#SBATCH -J genotype_filter_vcfs
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

# a pipeline to map trimmed fastq reads to a reference genome, produce cleaned, sorted, quality filtered bam files, and call snps using both gatk HaplotypeCaller and samtools mpileup
# a bash and conda-dependent pipeline

# usage: sbatch scripts/genotype_vcfs.sbatch

# define global variables
wd="/project/noujdine_61/jchancel/gallo_oa_popgen_pipeline/"
outdir="${wd}/var-calling/alignment/gatk"
samples_file="${wd}/bwa_samples.txt"
sample_id="$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)"
index="/project/noujdine_61/jchancel/indices/mussel/mytilus/galloprovinciales/MGAL_10.fa"
dictionary="/project/noujdine_61/jchancel/indices/mussel/mytilus/galloprovinciales/MGAL_10.dict"
rawvcf="${outdir}/gallo_OA.vcf.gz"
rawsnps="${outdir}/raw_snps_recal.vcf.gz"
snps="${outdir}/gallo_snps_passed_sites.vcf.gz"
biallelicsnps="${outdir}/biallelic_gallo_snps_passed_sites.recode.vcf"

# activate conda module and environments

module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
module load gcc/11.3.0
module load picard/2.26.2
module load gatk
module load samtools

conda activate gatk

# create accompanying index files for gvcfs

gatk IndexFeatureFile -I ${outdir}/MG22_D2_1.g.vcf.gz
gatk IndexFeatureFile -I ${outdir}/MG22_D2_2.g.vcf.gz
gatk IndexFeatureFile -I ${outdir}/MG22_D2_3.g.vcf.gz
gatk IndexFeatureFile -I ${outdir}/MG22_D23_C_1.g.vcf.gz
gatk IndexFeatureFile -I ${outdir}/MG22_D23_C_2.g.vcf.gz
gatk IndexFeatureFile -I ${outdir}/MG22_D23_C_3.g.vcf.gz
gatk IndexFeatureFile -I ${outdir}/MG22_D25_T_1.g.vcf.gz
gatk IndexFeatureFile -I ${outdir}/MG22_D25_T_2.g.vcf.gz
gatk IndexFeatureFile -I ${outdir}/MG22_D25_T_3.g.vcf.gz

# combine gvcfs from haplotype caller into single gvcf file with gatk

echo "Combining gvcfs from haplotype caller into single gvcf file with gatk."

gatk --java-options "-Xmx150g" CombineGVCFs \
-R ${index} \
--sequence-dictionary ${dictionary} \
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

echo "Genotyping ${outdir}/gallo_cohort.g.vcf.gz with GATK to obtain single, cohort vcf file."

gatk --java-options "-Xmx150g" GenotypeGVCFs \
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

# Extract only filtered, passed sites into hard-filtered, whole-cohort vcf file

echo "Hard filtering ${outdir}/filtered_snps.vcf.gz with GATK."

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

conda deactivate

# filter for only biallelic sites using vcftools
module purge
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate vcftools
module load gcc/11.3.0
module load gcc/8.3.0
module load gcc/9.2.0
module load intel/18.0.4
module load intel/19.0.4
module load vcftools/0.1.14

echo "Filtering ${snps} for biallelic sites with vcftools."

vcftools --gzvcf ${snps} --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out ${outdir}/biallelic_gallo_snps_passed_sites

conda deactivate

# output sample-specific vcfs

echo "Outputting sample-specific vcf files."

module purge
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
module load gcc/11.3.0
module load picard/2.26.2
module load gatk
module load samtools

conda activate gatk

gatk SelectVariants \
-V ${biallelicsnps} \
-sn MG22_D23_C_1 \
-sn MG22_D23_C_2 \
-sn MG22_D23_C_3 \
-O ${outdir}/MG22_D23_C_biallelic_snps_passedfilteredsites.vcf.gz

gatk SelectVariants \
-V ${biallelicsnps} \
-sn MG22_D25_T_1 \
-sn MG22_D25_T_2 \
-sn MG22_D25_T_3 \
-O ${outdir}/MG22_D25_T_biallelic_snps_passedfilteredsites.vcf.gz

gatk SelectVariants \
-V ${biallelicsnps} \
-sn MG22_D2_1 \
-sn MG22_D2_2 \
-sn MG22_D2_3 \
-O ${outdir}/MG22_D2_biallelic_snps_passedfilteredsites.vcf.gz

conda deactivate

echo "Disco!"