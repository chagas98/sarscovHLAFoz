#!/bin/bash -ue
mkdir -p B.1.1.28

cd B.1.1.28

# Run bcftools to perform variant calling on the input file using the reference file

bcftools mpileup --max-depth 3000 -f ../NC_045512.2.fasta ../B.1.1.28.bam | bcftools call -c --ploidy 1 -Ov -o B.1.1.28.vcf
