#!/bin/bash -ue
mkdir -p BE.1

cd BE.1

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../BE.1.fasta --secondary=no -N 0 > BE.1.sam

# Sort the aligned SAM file and save as BAM file
samtools sort BE.1.sam -o BE.1.bam

# Print the sorting status of the BAM file
samtools stats BE.1.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index BE.1.bam  BE.1.bam.bai
