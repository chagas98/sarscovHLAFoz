#!/bin/bash -ue
mkdir -p BE.9

cd BE.9

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../BE.9.fasta --secondary=no -N 0 > BE.9.sam

# Sort the aligned SAM file and save as BAM file
samtools sort BE.9.sam -o BE.9.bam

# Print the sorting status of the BAM file
samtools stats BE.9.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index BE.9.bam  BE.9.bam.bai
