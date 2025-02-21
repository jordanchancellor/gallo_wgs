#find -name "*.fq.gz" -execdir echo {} ';' | cut -c 3- | sed 's/\..*//g' | sed 's/_[12]//g' | sort -u > samples_file.txt
#ls */*/*.fq.gz | sed 's/.*\///g' | sed 's/\.fq\.gz//g' | sed 's/_[12]//g' | sort -u > samples_file.txt
ls 01.RawData/*/*.fq.gz | sed 's/.*\///g' | sed 's/\.fq\.gz//g' | rev | cut -c 3- | rev | sort -u > samples_file.txt
