#!/bin/bash -ue
mkdir -p BA.5.3

cd BA.5.3

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../BA.5.3.fasta --secondary=no -N 0 > BA.5.3.sam

# Sort the aligned SAM file and save as BAM file
samtools sort BA.5.3.sam -o BA.5.3.bam

# Print the sorting status of the BAM file
samtools stats BA.5.3.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index BA.5.3.bam  BA.5.3.bam.bai
