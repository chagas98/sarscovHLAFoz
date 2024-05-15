#!/bin/bash -ue
mkdir -p BE.10

cd BE.10

# Run bcftools to perform variant calling on the input file using the reference file

bcftools mpileup --max-depth 3000 -f ../NC_045512.2.fasta ../BE.10.bam | bcftools call -c --ploidy 1 -Ov -o BE.10.vcf
