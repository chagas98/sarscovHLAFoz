#!/bin/bash -ue
mkdir -p P.2

cd P.2

# Run bcftools to perform variant calling on the input file using the reference file

bcftools mpileup --max-depth 3000 -f ../NC_045512.2.fasta ../P.2.bam | bcftools call -c --ploidy 1 -Ov -o P.2.vcf
