#!/bin/bash -ue
mkdir -p B.1.1.28

cd B.1.1.28

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../B.1.1.28.fasta --secondary=no -N 0 > B.1.1.28.sam

# Sort the aligned SAM file and save as BAM file
samtools sort B.1.1.28.sam -o B.1.1.28.bam

# Print the sorting status of the BAM file
samtools stats B.1.1.28.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index B.1.1.28.bam  B.1.1.28.bam.bai
