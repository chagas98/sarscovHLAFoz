#!/bin/bash -ue
mkdir -p DL.1

cd DL.1

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../DL.1.fasta --secondary=no -N 0 > DL.1.sam

# Sort the aligned SAM file and save as BAM file
samtools sort DL.1.sam -o DL.1.bam

# Print the sorting status of the BAM file
samtools stats DL.1.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index DL.1.bam  DL.1.bam.bai
