#!/bin/bash
#SBATCH --mem=50GB
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=12
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH -J mpileup
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

# fix mate information on bam files with picard
# bash conda-dependent

# usage: sbatch --array=1-9 scripts/samtools_mpileup.sbatch

# define global variables and checkpoints

wd="/project/noujdine_61/jchancel/gallo_oa_popgen_pipeline/"
outdir="${wd}/var-calling/alignment/picard/"
popooldir="${wd}/var-calling/popoolation/"
index="/project/noujdine_61/jchancel/indices/mussel/mytilus/galloprovinciales/MGAL_10.fa"
popoolsoftware="/project/noujdine_61/jchancel/software/popoolation2_1201/"
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
# samfile="${wd}/var-calling/alignment/${sample_id}.sam"
# bamfile="${outdir}/${sample_id}.bam"
# sortedbamfile="${outdir}/${sample_id}.sorted.bam"
# dupesremovedbamfile="${outdir}/${sample_id}_duplicates_removed.sorted.bam"
# dupesremovedqfbamfile="${outdir}/${sample_id}_duplicates_removed.qf.sorted.bam"
# fixedmatesbam="${outdir}/${sample_id}_duplicates_removed.qf.matesfixed.sorted.bam"
cleanbam="${outdir}/${sample_id}_duplicates_removed.qf.matesfixed.sorted.calmd.bam"
mpileupfile="${popooldir}/${sample_id}.mpileup"
mpileupsync="${popooldir}/${sample_id}.mpileup.sync"
gtffile="${popooldir}/${sample_id}.mpileup.sync.indels.gtf"


# activate conda module and environments
module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
conda activate samtools

# compile snps by creating samtools mpileup file

echo "Creating mpileup file from ${cleanbam} with samtools."

samtools mpileup -B -Q 0 -f ${index} ${cleanbam} > ${popooldir}/${sample_id}.mpileup

conda deactivate

# synchronize file with popoolation 

echo "Syncing ${mpileupfile} with popoolation."

perl ${popoolsoftware}/mpileup2sync.pl --fastq-type illumina --min-qual 20 --input ${mpileupfile} --output ${popooldir}/${sample_id}.mpileup.sync

# filter indels with popoolation

echo "filtering indels on ${mpileupsync} with popoolation."

perl ${popoolsoftware}/indel_filtering/identify-indel-regions.pl \
--input ${mpileupsync} \
--output ${popooldir}/${sample_id}.mpileup.sync.indels.gtf \
--indel-window 5

# filter by gtf with popoolation

echo "filtering indels by gtf on ${gtffile} with popoolation."

perl ${popoolsoftware}/indel_filtering/filter-sync-by-gtf.pl \
--gtf ${gtffile} \
--input ${mpileupsync} \
--output ${popooldir}/${sample_id}_idf.mpileup.sync

# subsample to uniform coverage with popoolation

echo "Subsampling ${mpileupsync} to univorm coverate with popoolation."

perl ${popoolsoftware}/subsample-synchronized.pl \
--input ${popooldir}/${sample_id}_idf.mpileup.sync \
--output ${popooldir}/${sample_id}_ss10_idf.mpileup.sync \
--target-coverage 10 \
--max-coverage '2%' \
--method withoutreplace

echo "Disco!"