#!/bin/bash -ue
mkdir -p BA.2.23

cd BA.2.23

# Run bcftools to perform variant calling on the input file using the reference file

bcftools mpileup --max-depth 3000 -f ../NC_045512.2.fasta ../BA.2.23.bam | bcftools call -c --ploidy 1 -Ov -o BA.2.23.vcf
