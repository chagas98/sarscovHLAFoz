#!/bin/bash

# Define the input directory and file names
input_directory="variants"
multiple_fasta="test.fasta"
config_file="config.txt"
reference_sequence="ref_sars.fasta"

# Create a directory to store modified FASTA files
mkdir -p $input_directory

# Colocar cada sequencia em um folder, remover brazil|, etc

awk -v input_dir="$input_directory" 'BEGIN{RS=">";FS="\n"} NR>1{
    fnme=$1".fasta"; 
    gsub(/Brazil\|/, "", fnme); 
    gsub(/\|/, "+", fnme); 
    foldername=gensub(/\.fasta/, "", "g", input_dir "/" fnme); 
    system("mkdir -p " foldername); 
    print ">" $0 > foldername "/" fnme; 
    close(fnme);
    }' $multiple_fasta

# Run minimap2, bcftools mpileup, and snpeff in each folder

cd $input_directory

for sequence_file in */*.fasta; do
    relative_path=$(dirname "$sequence_file")
    sequence_name=$(basename "$sequence_file")

    #Input Fasta File
    query_fasta_file="$relative_path/$sequence_name"
    
    #Output Files
    bam_file="$( echo ${query_fasta_file} | sed 's/.fasta/.bam/')"
    vcf_file="$( echo ${query_fasta_file} | sed 's/.fasta/.vcf/')"
    snpeff_file="$( echo ${query_fasta_file} | sed 's/.fasta/.snpeff.vcf/')"
    
    echo $query_fasta_file
    
    # Run minimap2
    minimap2 -a -t 8 ../$reference_sequence $query_fasta_file --secondary=no > $bam_file

    # Run bcftools mpileup and call variants
    bcftools mpileup --max-depth 3000 -f ../$reference_sequence $bam_file | bcftools call -c --ploidy 1 -Ov -o $vcf_file

    # Run snpeff
    snpEff ann -v -c ../snpEff.config -csvStats ${query_fasta_file}_stats.csv -s ${query_fasta_file}_stats.html NC_045512.2 $vcf_file > $snpeff_file    

done

cd ..

