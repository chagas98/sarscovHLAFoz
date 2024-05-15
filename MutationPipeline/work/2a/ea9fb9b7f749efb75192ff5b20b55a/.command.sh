#!/bin/bash -ue
mkdir -p BA.4.6

cd BA.4.6

# Run bcftools to perform variant calling on the input file using the reference file

bcftools mpileup --max-depth 3000 -f ../NC_045512.2.fasta ../BA.4.6.bam | bcftools call -c --ploidy 1 -Ov -o BA.4.6.vcf
