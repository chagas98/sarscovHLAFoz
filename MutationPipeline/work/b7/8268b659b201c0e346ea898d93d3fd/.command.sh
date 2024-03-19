#!/bin/bash -ue
mkdir -p AY.122

cd AY.122

# Run Minimap2 to align the input file to the reference file
minimap2 -a ../NC_045512.2.fasta ../AY.122.fasta --secondary=no -N 0 > AY.122.sam

# Sort the aligned SAM file and save as BAM file
samtools sort AY.122.sam -o AY.122.bam

# Print the sorting status of the BAM file
samtools stats AY.122.bam | grep "is sorted:"

# Create an index file for the sorted BAM file
samtools index AY.122.bam  AY.122.bam.bai
