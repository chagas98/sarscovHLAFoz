#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Defines the pipeline input parameters (with a default value for each one).
 * Each of the following parameters can be specified as command line options.
 */
params.input_fasta = "$baseDir/test/test.fasta"
params.out = "$baseDir/result.txt"
params.reference = "$baseDir/test/ref_sars.fasta"
params.dir_results = "$baseDir/results"
params.n_perc = 5
params.dir_mutations = "mutations"
params.snpeff_db = 'NC_045512.2'
params.snpEff_config = "$baseDir/config/snpEff.config"
params.renviron = "$baseDir/config/.Renviron"


// Split Fasta into lineages and filter quality process
process getGISAIDseq {

    publishDir "$params.dir_results/raw_data", mode: 'copy', pattern: '**.{fasta,csv}', overwrite: true

    conda '/home/samuel/anaconda3/envs/TCC'

    input:
    path renviron

    output:
    path '**.fasta', emit: gisaid_fasta
    path '**.csv', emit: gisaid_csv

    script:
    """
    getvariants.R $renviron
    """
}

// Split Fasta into lineages and filter quality process
process clusteringByLineage {

    publishDir "$params.dir_results/$params.dir_mutations", mode: 'copy', pattern: '**.fasta', overwrite: true

    conda '/home/samuel/anaconda3/envs/TCC'

    input:
    path input
    val quality

    output:
    path '**.fasta'

    script:
    """
    filtering.py --input_fasta $input --perc_n $quality
    """
}

// Define the Minimap2 process
process minimap2 {

    publishDir "$params.dir_results/$params.dir_mutations/${input.baseName}", mode: 'copy', pattern:  '**.{bam,sam}', overwrite: true

    // Set Conda env
    conda '/home/samuel/anaconda3/envs/TCC'

    input:
    path input
    path reference

    output:
    path "**.sam", emit: aln_sam
    path "**.bam", emit: aln_sorted_bam

    script:
    """
    minimap2 -a $reference $input --secondary=no -N 0 > ${input.baseName}.sam
    
    samtools sort ${input.baseName}.sam -o ${input.baseName}.bam

    samtools stats ${input.baseName}.bam | grep "is sorted:"
    """
}

// Define the Minimap2 process
process variantcaller {

    publishDir "$params.dir_results/$params.dir_mutations/${input.baseName}", mode: 'copy', pattern: "**.vcf", overwrite: true

    // Set Conda env
    conda '/home/samuel/anaconda3/envs/TCC'

    input:
    path input
    path reference

    output:
    path "**.vcf" 

    script:
    """
    bcftools mpileup --max-depth 3000 -f $reference $input | bcftools call -c --ploidy 1 -Ov -o ${input.baseName}.vcf
    """
}

// Define the Minimap2 process
process annotationSnpeff {

    publishDir "$params.dir_results/$params.dir_mutations/${input.baseName}", mode: 'copy', pattern: '**.{vcf,csv,html}', overwrite: true

    // Set Conda env
    conda '/home/samuel/anaconda3/envs/TCC'

    input:
    path input
    val snpeff_db

    output:
    path "**.vcf", emit: vcf
    path "**.csv", emit: csv
    path "**.html", emit: html

    script:
    """
    snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html $snpeff_db $input > ${input.baseName}.snpeff.vcf
    """
}


workflow {
    /*
     * Create a channel emitting the given query fasta file(s).
     * Split the file into chunks containing as many sequences as defined by the parameter 'chunkSize'.
     * Finally, assign the resulting channel to the variable 'ch_fasta'
     */

    getGISAIDseq(params.renviron)

    split_fasta_ch = clusteringByLineage(getGISAIDseq.out.gisaid_fasta, params.n_perc)

    minimap2(split_fasta_ch.flatten(), params.reference)

    bcftools_ch = variantcaller(minimap2.out.aln_sorted_bam, params.reference)

    annotationSnpeff(bcftools_ch, params.snpeff_db)

}

