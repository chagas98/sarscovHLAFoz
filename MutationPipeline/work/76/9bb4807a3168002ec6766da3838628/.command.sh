#!/bin/bash -ue
mkdir -p BA.1.5

cd BA.1.5

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../BA.1.5.fasta --secondary=no -N 0 > BA.1.5.sam

# Sort the aligned SAM file and save as BAM file
samtools sort BA.1.5.sam -o BA.1.5.bam

# Print the sorting status of the BAM file
samtools stats BA.1.5.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index BA.1.5.bam  BA.1.5.bam.bai