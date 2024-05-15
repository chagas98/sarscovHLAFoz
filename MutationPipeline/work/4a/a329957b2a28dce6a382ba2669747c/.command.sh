#!/bin/bash -ue
mkdir -p BE.9

cd BE.9

# Run bcftools to perform variant calling on the input file using the reference file

bcftools mpileup --max-depth 3000 -f ../NC_045512.2.fasta ../BE.9.bam | bcftools call -c --ploidy 1 -Ov -o BE.9.vcf
