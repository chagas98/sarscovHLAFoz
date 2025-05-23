#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Define the default parameters
 */ 
params.n_perc = 5
params.n_seqs = 1
params.ref_seq = 'NC_045512.2'
params.snpEff_config = "$baseDir/config/snpEff.config"
params.renviron = "$baseDir/config/.Renviron"
params.versions = "$baseDir/config/versions.csv"
params.start_date = "2020-01-01"
params.end_date = "2023-01-01"
params.city_name = "'South America / Brazil / Parana / Foz'"
params.max_length = 14
params.min_length = 10
params.alleles = "B*07:02;B*44:03;B*15:01;B*51:01;B*18:01;B*35:01;B*38:01;B*08:01;B*49:01;B*14:02;B*39:01;B*27:05;B*40:01;B*53:01;B*57:01"

params.predict_tool = "netmhcpan-4.1"

/* 
 * Import modules 
 */
include {
  GET_REFSEQ_SARS;
  GET_GISAID_SEQUENCES;
  SELECT_BY_LINEAGE;
  RUN_MAPPING_MINIMAP2;
  VARIANT_CALLER;
  RUN_ANNOT_SNPEFF;
  RUN_EPYTOPE_PREDICTION;
  PLOTTING} from './modules.nf' 

/* 
 * main pipeline logic
 */
workflow {

    // Get SARS-CoV-2 RefSeq
    GET_REFSEQ_SARS(params.ref_seq)

    // Get GISAID sequences
    GET_GISAID_SEQUENCES(params.renviron, 
                         params.start_date, 
                         params.end_date, 
                         params.city_name)

    // Select sequences by lineage
    SELECT_BY_LINEAGE(GET_GISAID_SEQUENCES.out.gisaid_fasta, 
                      params.n_perc, 
                      params.n_seqs)

    // Run mapping with minimap2
    RUN_MAPPING_MINIMAP2(SELECT_BY_LINEAGE.out.fasta.flatten(), 
                         GET_REFSEQ_SARS.out)

    // Perform variant calling
    VARIANT_CALLER(RUN_MAPPING_MINIMAP2.out.aln_sorted_bam, 
                   GET_REFSEQ_SARS.out)

    // Run annotation with snpEff
    RUN_ANNOT_SNPEFF(VARIANT_CALLER.out, 
                     params.ref_seq)

    // Run epytope prediction
    RUN_EPYTOPE_PREDICTION(RUN_ANNOT_SNPEFF.out.vcf,
                           params.versions,
                           params.max_length,
                           params.min_length,
                           params.alleles,
                           params.predict_tool)

}

