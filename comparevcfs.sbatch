#!/bin/bash
#SBATCH --mem=0
#SBATCH --time=1-0
#SBATCH --partition=main
#SBATCH --cpus-per-task=12
#SBATCH --error=%x.err
#SBATCH --output=%x.out
#SBATCH -J comparevcfs
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

# define global variables
wd="/project/noujdine_61/jchancel/gallo_oa_popgen_pipeline/"
outdir="${wd}/var-calling/alignment/gatk"
day2vcf="${outdir}/MG22_D23_C_biallelic_snps_passedfilteredsites.vcf.gz"
controlvcf="${outdir}/MG22_D25_T_biallelic_snps_passedfilteredsites.vcf.gz"
treatmentvcf="${outdir}/MG22_D2_biallelic_snps_passedfilteredsites.vcf.gz"

# activate conda environment and load required modules

module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
module load gcc/11.3.0
module load intel/19.0.4
module load openblas/0.3.21
module load vcftools/0.1.14
module load bcftools/1.14

conda activate vcftools

# compare vcf files for differing sites

echo "Comparing differing variant sites between treatment and control cohorts with vcftools."

vcftools --gzvcf ${controlvcf} --gzdiff ${treatmentvcf} --diff-indv-discordance --out D23_C_v_D25_T_diff_results

conda deactivate

conda activate bcftools

echo "Extracting variant sites unique to each treatment with bcftools."

bcftools isec -p D23_C_v_D25_T_isec_output -n-1 -c snps ${controlvcf} ${treatmentvcf}

echo "Disco!"
