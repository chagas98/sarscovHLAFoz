#!/bin/bash -ue
mkdir -p AY.46.3

cd AY.46.3

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../AY.46.3.fasta --secondary=no -N 0 > AY.46.3.sam

# Sort the aligned SAM file and save as BAM file
samtools sort AY.46.3.sam -o AY.46.3.bam

# Print the sorting status of the BAM file
samtools stats AY.46.3.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index AY.46.3.bam  AY.46.3.bam.bai
