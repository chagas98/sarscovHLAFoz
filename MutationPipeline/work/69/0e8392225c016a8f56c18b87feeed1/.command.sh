#!/bin/bash -ue
variant=B.1.1

mkdir -p ${variant}

cd ${variant}

#mkdir -p prediction

#cd prediction

export PYTHONWARNINGS="ignore::RuntimeWarning"
export PYTHONWARNINGS="ignore::FutureWarning"

# Run epaa_mod.py to predict epitope for the input VCF file    
epaa_mod.py --identifier ${variant} --alleles "B*07:02;B*44:03;B*15:01;B*51:01;B*18:01;B*35:01;B*38:01;B*08:01;B*49:01;B*14:02;B*39:01;B*27:05;B*40:01;B*53:01;B*57:01" --tools "netmhcpan-4.1" --max_length 14 --min_length 10 --versions ../versions.csv --variant_lineage ${variant} --somatic_mutation ../B.1.1.snpEff.vcf
