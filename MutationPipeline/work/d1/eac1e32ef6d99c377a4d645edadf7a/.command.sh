#!/bin/bash -ue
mkdir -p BQ.1.1

cd BQ.1.1

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../BQ.1.1.fasta --secondary=no -N 0 > BQ.1.1.sam

# Sort the aligned SAM file and save as BAM file
samtools sort BQ.1.1.sam -o BQ.1.1.bam

# Print the sorting status of the BAM file
samtools stats BQ.1.1.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index BQ.1.1.bam  BQ.1.1.bam.bai
