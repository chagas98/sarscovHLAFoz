#!/bin/bash -ue
mkdir -p BQ.1.1

cd BQ.1.1

# Run bcftools to perform variant calling on the input file using the reference file

bcftools mpileup --max-depth 3000 -f ../NC_045512.2.fasta ../BQ.1.1.bam | bcftools call -c --ploidy 1 -Ov -o BQ.1.1.vcf
