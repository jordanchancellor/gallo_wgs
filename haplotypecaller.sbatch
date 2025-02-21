#!/bin/bash
#SBATCH --mem=50g
#SBATCH --partition=oneweek
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=12
#SBATCH --error=outfiles/gatk2/%x_%a.err
#SBATCH --output=outfiles/gatk2/%x_%a.out
#SBATCH -J haplotypecaller
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

index="/scratch2/jchancel/gallo_oa_popgen/var-calling/reference_data/GCA_025277285.1_MgalMED_genomic.fna"
outdir="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/gatk2/"
samples_file="/scratch2/jchancel/gallo_oa_popgen/bwa_samples.txt"
sample_id="$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)"
samfile="var-calling/alignment/bwa/${sample_id}.sam"
bamfile="var-calling/alignment/bwa/${sample_id}.sorted.bam"
cleanbam="var-calling/alignment/gatk2/${sample_id}.cleaned.sorted.bam"
dedupbamfile="var-calling/alignment/gatk2/${sample_id}_marked_duplicates.bam"
targetint="var-calling/alignment/gatk/${sample_id}.realignertargetcreator.intervals"
dict="var-calling/reference_data/GCA_025277285.1_MgalMED_genomic.dict"

#activate conda environment and load required modules
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate samtools
module load gcc/11.3.0
module load samtools
module load picard/2.26.2
module load gatk

#index fixed_mate bam files

echo "Indexing ${outdir}/${sample_id}.fixed_mate.bam using samtools."
samtools index ${outdir}/${sample_id}.fixed_mate.bam

conda deactivate
conda activate gatk

echo "Calling variants using Haplotype caller."

gatk  HaplotypeCaller \
-R ${index} \
-I ${outdir}/${sample_id}.fixed_mate.bam \
-O ${outdir}/${sample_id}.g.vcf.gz \
-ERC GVCF
