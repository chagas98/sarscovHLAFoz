#!/bin/bash -ue
mkdir -p AY.101

cd AY.101

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../AY.101.fasta --secondary=no -N 0 > AY.101.sam

# Sort the aligned SAM file and save as BAM file
samtools sort AY.101.sam -o AY.101.bam

# Print the sorting status of the BAM file
samtools stats AY.101.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index AY.101.bam  AY.101.bam.bai
