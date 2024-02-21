#!/bin/bash -ue
mkdir -p AY.99.2

cd AY.99.2

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../AY.99.2.fasta --secondary=no -N 0 > AY.99.2.sam

# Sort the aligned SAM file and save as BAM file
samtools sort AY.99.2.sam -o AY.99.2.bam

# Print the sorting status of the BAM file
samtools stats AY.99.2.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index AY.99.2.bam  AY.99.2.bam.bai
