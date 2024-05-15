#!/bin/bash -ue
mkdir -p BA.2.23

cd BA.2.23

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../BA.2.23.fasta --secondary=no -N 0 > BA.2.23.sam

# Sort the aligned SAM file and save as BAM file
samtools sort BA.2.23.sam -o BA.2.23.bam

# Print the sorting status of the BAM file
samtools stats BA.2.23.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index BA.2.23.bam  BA.2.23.bam.bai
