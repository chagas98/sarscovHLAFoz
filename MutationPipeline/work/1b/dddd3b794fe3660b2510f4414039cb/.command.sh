#!/bin/bash -ue
mkdir -p DL.1

cd DL.1

# Run bcftools to perform variant calling on the input file using the reference file

bcftools mpileup --max-depth 3000 -f ../NC_045512.2.fasta ../DL.1.bam | bcftools call -c --ploidy 1 -Ov -o DL.1.vcf
