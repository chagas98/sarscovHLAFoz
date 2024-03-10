#!/bin/bash -ue
mkdir -p P.7

cd P.7

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../P.7.fasta --secondary=no -N 0 > P.7.sam

# Sort the aligned SAM file and save as BAM file
samtools sort P.7.sam -o P.7.bam

# Print the sorting status of the BAM file
samtools stats P.7.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index P.7.bam  P.7.bam.bai
