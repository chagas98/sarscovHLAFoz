## Install Dependencies

- Create a conda environment

conda create -n TCC

- Activate the conda environment

conda activate TCC

- Install minimap2

conda install -c bioconda minimap2

- Install bcftools

conda install -c bioconda bcftools

- Install vcftools

conda install -c bioconda vcftools

- Install snpeff

conda install -c bioconda snpeff


## Find Mutations

    minimap2 -ax asm10 ref_sars.fasta test.fasta > aln.bam

    samtools stats aln.bam | grep "is sorted:"

    samtools sort aln.bam -o aln_sorted.bam

    bcftools mpileup --max-depth 3000 -f ref_sars.fasta aln_sorted.bam | bcftools call -c --ploidy 1 -Ov -o call_variants.vcf

    cp /home/samuel/anaconda3/envs/TCC/share/snpeff-5.2-0/snpEff.config snpEff.config

    snpEff ann -v -c snpEff.config -s annotated_vcf.vcf NC_045512.2 call_variants.vcf > annotated_vcf_ann.vcf

    snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv \-s annotated_vcf_ann.html NC_045512.2 call_variants.vcf > annotated_vcf_ann.vcf











