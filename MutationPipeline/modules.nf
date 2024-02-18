

// Get RefSeq from NCBI
process GET_REFSEQ_SARS {

    input:
    val seqId

    output:
    path '**.fasta'

    script:
    """
    # Download the fasta file from NCBI using the provided sequence ID

    wget -O ${seqId}.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&report=fasta&id=${seqId}"
    """
}

// Get GISAID sequences
process GET_GISAID_SEQUENCES {

    input:
    path renviron
    val start_date
    val end_date
    val city_name

    output:
    path '**.fasta', emit: gisaid_fasta
    path '**.csv', emit: gisaid_csv

    script:
    """
    # Run the getvariants.R script using the provided renviron file, start-to-end date and city name
    
    getvariants.R $renviron $start_date $end_date $city_name
    """
}

// Split Fasta into lineages and filter quality process
process SELECT_BY_LINEAGE {

    input:
    path input
    val quality
    val number_seqs

    output:
    path '**.fasta', emit: fasta
    path '**.csv', emit: csv

    script:
    """
    # Filter the input fasta file based on quality and number of sequences

    filtering.py --input_fasta $input --perc_n $quality --minimum $number_seqs
    """
}

// Define the Minimap2 process
process RUN_MAPPING_MINIMAP2 {

    input:
    path input
    path reference

    output:
    path "**.sam", emit: aln_sam
    path "**.bam", emit: aln_sorted_bam
    path "**.bai", emit: aln_sorted_index_bam

    script:
    """
    mkdir -p ${input.baseName}
    
    cd ${input.baseName}

    # Run Minimap2 to align the input file to the reference file
    minimap2 -a ../$reference ../$input --secondary=no -N 0 > ${input.baseName}.sam
    
    # Sort the aligned SAM file and save as BAM file
    samtools sort ${input.baseName}.sam -o ${input.baseName}.bam

    # Print the sorting status of the BAM file
    samtools stats ${input.baseName}.bam | grep "is sorted:"

    # Create an index file for the sorted BAM file
    samtools index ${input.baseName}.bam  ${input.baseName}.bam.bai
    """
}

// Variants call
process VARIANT_CALLER {

    input:
    path input
    path reference

    output:
    path "**.vcf" 

    script:
    """
    mkdir -p ${input.baseName}
    
    cd ${input.baseName}

    # Run bcftools to perform variant calling on the input file using the reference file

    bcftools mpileup --max-depth 3000 -f ../$reference ../$input | bcftools call -c --ploidy 1 -Ov -o ${input.baseName}.vcf
    """
}

// Define the annotation of the variants
process RUN_ANNOT_SNPEFF {

    input:
    path input
    val snpeff_db

    output:
    path "**.vcf", emit: vcf
    path "**.csv", emit: csv
    path "**.html", emit: html

    script:
    """

    mkdir -p ${input.baseName}
    
    cd ${input.baseName}

    # Run snpEff to annotate the input VCF file using the provided snpeff_db
    snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html $snpeff_db ../$input > ${input.baseName}.snpeff.vcf
    """
}


process RUN_EPYTOPE_PREDICTION {

    input:
    path folder
    path input
    path versions

    output:
    path "**.csv", emit: csv

    script:
    """
    
    cd ${folder.baseName}

    mkdir -p prediction

    cd prediction

    # Run epaa_mod.py to predict epitope for the input VCF file    
    epaa_mod.py --identifier ${input.baseName} --alleles 'A*01:01;A*02:01' --tools 'netmhcpan-4.1' --max_length 15 --min_length 12 --versions $versions --variant_lineage ${input.baseName} --somatic_mutation $input
    """
}
