#!/bin/bash
#SBATCH --mem=50gb
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=12
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH -J trimmingstats
#SBATCH --mail-type=END
#SBATCH --mail-user=jchancel@usc.edu

# usuage: sbatch scripts/trimmingstats.sbatch

wd="/project/noujdine_61/jchancel/gallo_oa_popgen_pipeline/"
dir="${wd}/gallo_trimmed/trimming_reports"
output_file="trimming_summary.csv"

# Define output CSV file and write the header
output_file="trimming_summary.csv"
echo "sample_id,total_reads_processed,reads_with_adapters,reads_written" > "$output_file"

# Loop through each input file as a command-line argument
for file in ${dir}/*.txt; do
    # Extract the sample ID by removing the suffix ".fq.gz_trimming_report.txt"
    sample_id=$(basename "$file" | sed 's/\.fq\.gz_trimming_report\.txt//')

    # Extract total reads processed
    total_reads=$(grep "Total reads processed:" "$file" | awk '{print $4}' | tr -d ',')

    # Extract reads with adapters
    reads_with_adapters=$(grep "Reads with adapters:" "$file" | awk '{print $4}' | tr -d ',')

    # Extract reads written (passing filters)
    reads_written=$(grep "Reads written (passing filters):" "$file" | awk '{print $5}' | tr -d ',')

    # Write the results to the output file
    echo "$sample_id,$total_reads,$reads_with_adapters,$reads_written" >> "$output_file"
done


