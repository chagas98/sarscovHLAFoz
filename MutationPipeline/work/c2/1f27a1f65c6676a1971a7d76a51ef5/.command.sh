#!/bin/bash -ue
mkdir -p P.2

cd P.2

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../P.2.fasta --secondary=no -N 0 > P.2.sam

# Sort the aligned SAM file and save as BAM file
samtools sort P.2.sam -o P.2.bam

# Print the sorting status of the BAM file
samtools stats P.2.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index P.2.bam  P.2.bam.bai
