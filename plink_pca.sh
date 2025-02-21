#!/bin/bash
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --partition=htcondor
#SBATCH --cpus-per-task=12
#SBATCH --error=outfiles/gatk2/%x.err
#SBATCH --output=outfiles/gatk2/%x.out
#SBATCH -J plinkpca
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

index="/scratch2/jchancel/gallo_oa_popgen/var-calling/reference_data/GCA_025277285.1_MgalMED_genomic.fna"
outdir="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/"
snps="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/gallo_snps_passed_sites.vcf.gz"

module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate plink

#make a BED file and the associated plink format files using the following command in plink:
#plink --const-fid 0 --vcf ${snps} --make-bed --allow-extra-chr --chr-set 14 no-xy --out ${outdir}/gallo_snps_allsamples

#carry out PCA using plink files
plink --pca --bfile ${outdir}/gallo_snps_allsamples --allow-extra-chr --chr-set 14 no-xy --out ${outdir}/gallo_snps_allsamples
