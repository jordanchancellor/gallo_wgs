#!/bin/bash
#SBATCH --mem=100g
#SBATCH --time=7-0
#SBATCH --partition=largemem
#SBATCH --nodelist=a02-10
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH -J popoolation
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

# a pipeline to map trimmed fastq reads to a reference genome, produce cleaned, sorted, quality filtered bam files, and call snps using both gatk HaplotypeCaller and samtools mpileup
# a bash and conda-dependent pipeline

# usage: sbatch scripts/popoolation.sbatch

# define global variables
wd="/project/noujdine_61/jchancel/gallo_oa_popgen_pipeline/"
mpileupoutdir="${wd}/var-calling/alignment/mpileup/"
popooldir="${wd}/var-calling/popoolation/"
popoolsoftware="/project/noujdine_61/jchancel/software/popoolation2_1201/"
index="/project/noujdine_61/jchancel/indices/mussel/mytilus/galloprovinciales/MGAL_10.fa"
mpileupfile="${mpileupoutdir}/control_treatment.mpileup"
mpileupsync="${popooldir}/control_treatment.mpileup.sync"

module load conda
source /spack/conda/miniconda3/23.3.1/etc/profile.d/conda.sh
module load gcc/11.3.0
module load samtools

# compile snps by creating samtools mpileup file

echo "Creating mpileup file for treatment and control cohorts with samtools."

conda activate samtools

samtools mpileup -B -Q 0 -f ${index} ${mpileupoutdir}/MG22_D23_C_1_duplicates_removed.sorted.bam ${mpileupoutdir}/MG22_D23_C_2_duplicates_removed.sorted.bam ${mpileupoutdir}/MG22_D23_C_3_duplicates_removed.sorted.bam ${mpileupoutdir}/MG22_D25_T_1_duplicates_removed.sorted.bam ${mpileupoutdir}/MG22_D25_T_2_duplicates_removed.sorted.bam ${mpileupoutdir}/MG22_D25_T_3_duplicates_removed.sorted.bam > ${mpileupoutdir}/control_treatment.mpileup

conda deactivate

# synchronize file with popoolation 

echo "Syncing ${mpileupfile} with popoolation."

perl ${popoolsoftware}/mpileup2sync.pl --fastq-type illumina --min-qual 20 --input ${mpileupfile} --output ${popooldir}/control_treatment.mpileup.sync

# calculate allele frequency differences

perl ${popoolsoftware}/snp-frequency-diff.pl --input ${popooldir}/control_treatment.mpileup.sync --output-prefix control_treatment --min-count 6 --min-coverage 50 --max-coverage 200

# calculate fst using sliding window approach

perl ${popoolsoftware}/fst-sliding.pl --input ${popooldir}/control_treatment.mpileup.sync --output control_treatment_w500.fst --min-count 6 --min-coverage 50 --max-coverage 200 --min-covered-fraction 1 --window-size 500 --step-size 500 --pool-size 500

# # index fst bam files

# conda activate samtools

# samtools index map/pop1.bam
# samtools index map/pop2.bam

# conda deactivate

# # convert file containing the Fst-values into the .igv format

# perl ${popoolsoftware}/export/pwc2igv.pl --input ${popooldir}/control_treatment.fst --output p1_p2.igv

# estimate the significance of allele frequency differences with Fisher's exact test

perl ${popoolsoftware}/fisher-test.pl --input ${popooldir}/control_treatment.mpileup.sync --output ${popooldir}/control_treatment.fet --min-count 6 --min-coverage 50 --max-coverage 200 --suppress-noninformative

# convert the Fisher's exact test results into the .igv file format:

perl ${popoolsoftware}/export/pwc2igv.pl --input ${popooldir}/control_treatment.fet --output ${popooldir}/control_treatment_fet.igv

#  detect consistent allele frequency changes in several biological replicates with Cochran-Mantel-Haenszel test

perl ${popoolsoftware}/cmh-test.pl --input ${popooldir}/control_treatment.mpileup.sync --output ${popooldir}/control_treatment.cmh --min-count 12 --min-coverage 50 --max-coverage 200 --population 1-4,2-5,3-6

# convert the results of the cmh-test into the .gwas format of IGV.

perl ${popoolsoftware}/export/cmh2gwas.pl --input ${popooldir}/control_treatment.cmh --output ${popooldir}/control_treatment.gwas --min-pvalue 1.0e-20