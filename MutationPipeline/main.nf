#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Defines the pipeline input parameters (with a default value for each one).
 * Each of the following parameters can be specified as command line options.
 */
params.query = "$baseDir/test/test.fasta"
params.out = "$baseDir/result.txt"
params.reference = "$baseDir/test/ref_sars.mmi"
params.dir= "$baseDir/variants"
params.chunkSize = 1
 
query_name = file(params.query).name
query_dir = file(params.query).parent
ref_fasta = file(params.reference).name

// Split Fasta into lineages and filter quality process
process clusteringByLineage {

    publishDir "variants/", mode: 'copy'

    conda '/home/samuel/anaconda3/envs/TCC'

    input:
    path query

    output:
    path 'output.log'

    script:
    """
    filtering.py --input_fasta $query
    """
}



// Define the Minimap2 process
process minimap2 {

    // Set Conda env
    conda '/home/samuel/anaconda3/envs/TCC'

    input:
    path query
    path reference

    output:
    path "aligned_${query.baseName}_sorted.bam" 

    script:
    """
    minimap2 -a $reference $query --secondary=no -N 0 > aligned_${query.baseName}.bam
    
    samtools sort aligned_${query.baseName}.bam -o aligned_${query.baseName}_sorted.bam

    samtools stats aligned_${query.baseName}_sorted.bam | grep "is sorted:"
    """
}

workflow {
    /*
     * Create a channel emitting the given query fasta file(s).
     * Split the file into chunks containing as many sequences as defined by the parameter 'chunkSize'.
     * Finally, assign the resulting channel to the variable 'ch_fasta'
     */

    Channel
        .fromPath(params.query)
        .splitFasta(by: params.chunkSize, file:true)
        .set { ch_fasta }

    split_fasta = clusteringByLineage(params.query)

    /*
     * Execute a Minimap2 alignment for each chunk emitted by the 'ch_fasta' channel
     * and emit the resulting Minimap2 alignments.
     */

    //ch_alignments = minimap2(ch_fasta, params.reference)

    /*
     * Collect all the sequences files into a single file
     * and print the resulting file contents when complete.
     */
    
}

