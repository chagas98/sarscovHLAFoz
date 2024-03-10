#!/bin/bash -ue
mkdir -p P.1.7

cd P.1.7

# Run bcftools to perform variant calling on the input file using the reference file

bcftools mpileup --max-depth 3000 -f ../NC_045512.2.fasta ../P.1.7.bam | bcftools call -c --ploidy 1 -Ov -o P.1.7.vcf
