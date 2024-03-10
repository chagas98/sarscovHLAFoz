#!/bin/bash -ue
mkdir -p P.1

cd P.1

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../P.1.fasta --secondary=no -N 0 > P.1.sam

# Sort the aligned SAM file and save as BAM file
samtools sort P.1.sam -o P.1.bam

# Print the sorting status of the BAM file
samtools stats P.1.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index P.1.bam  P.1.bam.bai
