#!/bin/bash -ue
mkdir -p BA.4

cd BA.4

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../BA.4.fasta --secondary=no -N 0 > BA.4.sam

# Sort the aligned SAM file and save as BAM file
samtools sort BA.4.sam -o BA.4.bam

# Print the sorting status of the BAM file
samtools stats BA.4.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index BA.4.bam  BA.4.bam.bai
