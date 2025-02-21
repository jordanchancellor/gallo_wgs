#!/bin/bash
#SBATCH --mem=0
#SBATCH --time=48:00:00
#SBATCH --partition=epyc-64
#SBATCH --cpus-per-task=12
#SBATCH --error=outfiles/gatk2/%x_%a.err
#SBATCH --output=outfiles/gatk2/%x_%a.out
#SBATCH -J extractandfilter
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

index="/scratch2/jchancel/gallo_oa_popgen/var-calling/reference_data/GCA_025277285.1_MgalMED_genomic.fna"
outdir="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/"
rawvcf="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/gallo_OA.vcf.gz"
snps="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/raw_snps_recal.vcf.gz"
indels="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/raw_indels_recal.vcf.gz"

#activate conda environment and load required modules
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate gatk
module load gcc/11.3.0
module load picard/2.26.2
module load gatk

#echo "Extracting SNPs and Indels"

gatk SelectVariants \
        -R ${index} \
        -V ${rawvcf} \
        -select-type-to-include SNP \
        -O ${outdir}/raw_snps_recal.vcf.gz

gatk SelectVariants \
        -R ${index} \
        -V ${rawvcf} \
        -select-type-to-include INDEL \
        -O ${outdir}/raw_indels_recal.vcf.gz

echo "Filtering SNPs"

gatk VariantFiltration \
-V ${snps} \
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
-O ${outdir}/filtered_snps_final.vcf.gz

echo "Filtering indels"

gatk VariantFiltration \
-V ${indels} \
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
-O ${outdir}/filtered_indels_final.vcf.gz

#echo "Exporting Variants to tab delimited table"

#gatk VariantsToTable \
#     -V scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk/filtered_snps_final.vcf \
#     -F CHROM -F POS -F TYPE -GF AD \
#     -O ${outdir}/filtered_snps_final.table

#gatk VariantsToTable \
#     -V scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk/filtered_indels_final.vcf \
#     -F CHROM -F POS -F TYPE -GF AD \
#     -O ${outdir}/filtered_indels_final.table

echo "Done."
