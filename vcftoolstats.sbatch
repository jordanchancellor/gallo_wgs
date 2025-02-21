#!/bin/bash
#SBATCH --mem=0
#SBATCH --time=48:00:00
#SBATCH --partition=htcondor
#SBATCH --cpus-per-task=12
#SBATCH --error=outfiles/gatk2/%x.err
#SBATCH --output=outfiles/gatk2/%x.out
#SBATCH -J vcftoolstats
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

index="/scratch2/jchancel/gallo_oa_popgen/var-calling/reference_data/GCA_025277285.1_MgalMED_genomic.fna"
outdir="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/"
snps="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/gallo_snps_passed_sites.vcf.gz"
controlvcf="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/MG22_D23_C_snps_passedfilteredsites.vcf.gz"
treatmentvcf="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/MG22_D25_T_snps_passedfilteredsites.vcf.gz"
day2vcf="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/MG22_D2_snps_passedfilteredsites.vcf.gz"
day2pop="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/day2popnames.txt"
controlpop="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/controlpopnames.txt"
treatmentpop="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/treatmentpopnames.txt"

#activate conda environment and load required modules
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate vcftools
module load gcc/11.3.0
module load gcc/8.3.0
module load gcc/9.2.0
module load intel/18.0.4
module load intel/19.0.4
module load vcftools/0.1.14

#calculate nucleotide diversity for subsetted vcfs
#vcftools --gzvcf ${day2vcf} --window-pi  10000 --out ${outdir}/MG22_D2
#vcftools --gzvcf ${treatmentvcf} --window-pi  10000 --out ${outdir}/MG22_D25_T
#vcftools --gzvcf ${controlvcf} --window-pi  10000 --out ${outdir}/MG22_D23_C

#calculate Tajima’s D
#vcftools --gzvcf ${day2vcf} --TajimaD 10000 --out ${outdir}/MG22_D2
#vcftools --gzvcf ${treatmentvcf} --TajimaD 10000 --out ${outdir}/MG22_D25_T
#vcftools --gzvcf ${controlvcf} --TajimaD 10000 --out ${outdir}/MG22_D25_T

#determine allele frequency for each population
#vcftools --gzvcf ${day2vcf} --freq --out ${outdir}/MG22_D2
#vcftools --gzvcf ${treatmentvcf} --freq --out ${outdir}/MG22_D25_T
#vcftools --gzvcf ${controlvcf} --freq --out ${outdir}/MG22_D25_T

#Estimate Fst population distance between control and treatment cohorts
#vcftools --gzvcf ${snps} --weir-fst-pop ${controlpop} --weir-fst-pop ${treatmentpop} --out ${outdir}/pop1_vs_2_FST 

#Estimage Fst population distance between day 2 and settlement days for both treatments
vcftools --gzvcf ${snps} --weir-fst-pop ${day2pop} --weir-fst-pop ${treatmentpop} --out ${outdir}/day2_vs_day25_T_FST
vcftools --gzvcf ${snps} --weir-fst-pop ${day2pop} --weir-fst-pop ${controlpop} --out ${outdir}/day2_vs_day23_C_FST
