#!/bin/bash
#SBATCH --mem=50g
#SBATCH --time=1-0
#SBATCH --partition=epyc-64
#SBATCH --cpus-per-task=12
#SBATCH --error=outfiles/bwa/%x_%a.err
#SBATCH --output=outfiles/bwa/%x_%a.out
#SBATCH -J bwa_sampe
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

index="/scratch2/jchancel/gallo_oa_popgen/var-calling/reference_data/GCA_025277285.1_MgalMED_genomic.fna"
trimmedreads="/scratch2/jchancel/gallo_oa_popgen/gallo_trimmed/"
outdir="/scratch2/jchancel/gallo_oa_popgen/var-calling/alignment/bwa2/"
samples_file="/scratch2/jchancel/gallo_oa_popgen/bwa_samples.txt"
sample_id="$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)"
reads1="${trimmedreads}/${sample_id}/${sample_id}_R1.fq.gz"
reads2="${trimmedreads}/${sample_id}/${sample_id}_R2.fq.gz"
sai1="${outdir}/${sample_id}_R1.sai"
sai2="${outdir}/${sample_id}_R2.sai"
header="$(zcat $reads1 | head -n 1)"
id="$(echo $header | awk '{print $1}' | sed 's/@//g' | cut -f 1-4 -d ":" | cut -f 3,4 -d ":" | sed 's/:/./g')"
sm="${sample_id}"
pu="${id}.${sm}"
lb="$(echo $header | awk '{print $2}' | cut -f 4 -d ":")"

echo "Read Group @RG\tID:$id\tSM:$sm\tLB:$lb\tPL:ILLUMINA"

#Load bwa conda environment
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate bwa

echo "Running bwa sampe on 2 read files associated with ${reads1} and ${reads2}."

cmd="bwa sampe -r $(echo "@RG\tID:$id\tSM:$sm\tLB:$lb\tPL:ILLUMINA") $index $sai1 $sai2 $reads1 $reads2"

echo "$cmd"

$cmd > "var-calling/alignment/bwa/${sample_id}.sam"

conda deactivate
conda activate samtools

echo "Converting sam files to sorted bam files with samtools."

samfile="var-calling/alignment/bwa/${sample_id}.sam"

cat ${samfile} | samtools view -bSh | samtools sort -o var-calling/alignment/bwa/${sample_id}.sorted.bam

bamfile="var-calling/alignment/bwa/${sample_id}.sorted.bam"

echo "Creating samtools index file from ${bamfile}"

samtools index -b ${bamfile} 

echo "Running idxstats on sorted, indexed bam file from ${bamfile}"

samtools idxstats ${bamfile} | tee ${outdir}/${sample_id}.idxstats.txt

echo "Creating mapping rate file with samtools flagstat from ${bamfile}"

samtools flagstat ${bamfile} > ${outdir}/${sample_id}.flagstat.txt

echo "Done."
