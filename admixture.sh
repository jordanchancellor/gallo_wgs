#!/bin/bash
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --partition=htcondor
#SBATCH --cpus-per-task=12
#SBATCH --error=outfiles/gatk2/%x.err
#SBATCH --output=outfiles/gatk2/%x.out
#SBATCH -J admixturek5
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

index="/scratch2/jchancel/gallo_oa_popgen/var-calling/reference_data/GCA_025277285.1_MgalMED_genomic.fna"
outdir="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/"
snps_numericChr="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/gallo_snps_passed_sites.numericChr.nocontigs.recode.vcf"
bedfile="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/gallo_snps_allsamples.numericChr.bed"

#load conda
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh

#make bed file for admixture using plink
conda activate plink
plink --const-fid 0 --vcf ${snps_numericChr} --geno 0.999 --make-bed --allow-extra-chr --chr-set 14 no-xy --out ${outdir}/gallo_snps_allsamples.numericChr
conda deactivate

#run admixture
conda activate admixture 
for K in {1..4}; do
admixture --cv ${bedfile} $K | tee log${K}.out; done

