#!/bin/bash -ue
mkdir -p BA.1.15

cd BA.1.15

# Run bcftools to perform variant calling on the input file using the reference file

bcftools mpileup --max-depth 3000 -f ../NC_045512.2.fasta ../BA.1.15.bam | bcftools call -c --ploidy 1 -Ov -o BA.1.15.vcf