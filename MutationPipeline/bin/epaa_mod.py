#!/usr/bin/env python
# Written by Christopher Mohr (adapted by Samyel Chagas) and released under the MIT license (2022).

import os
import sys
import logging
import csv
import re
import vcf
import argparse
import urllib
import itertools
import pandas as pd
import numpy as np
import epytope.Core.Generator as generator
import math
import json

from epytope.IO.MartsAdapter import MartsAdapter
from epytope.Core.Variant import Variant, VariationType, MutationSyntax
from epytope.EpitopePrediction import EpitopePredictorFactory
from epytope.IO.ADBAdapter import EIdentifierTypes
from epytope.IO.UniProtAdapter import UniProtDB
from epytope.Core.Allele import Allele
from epytope.Core.Peptide import Peptide
from epytope.Core.Protein import Protein
from epytope.Core.Transcript import Transcript
from Bio import SeqUtils
from datetime import datetime

__author__ = "Christopher Mohr (from Epytope python package) Adapted by Samuel Chagas de Assis"
VERSION = "1.1"

# instantiate global logger object
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)

#ID_SYSTEM_USED = EIdentifierTypes.ENSEMBL
ID_SYSTEM_USED = EIdentifierTypes.REFSEQ
transcriptProteinMap = {}
transcriptSwissProtMap = {}

_allowed_aas = frozenset('ACDEFGHIKLMNPQRSTVWY')

def get_epytope_annotation(vt, p, r, alt):
    if vt == VariationType.SNP:
        return p, r, alt
    elif vt == VariationType.DEL or vt == VariationType.FSDEL:
        # more than one observed ?
        if alt != "-":
            alternative = "-"
            reference = r[len(alt) :]
            position = p + len(alt)
        else:
            return p, r, alt
    elif vt == VariationType.INS or vt == VariationType.FSINS:
        if r != "-":
            position = p
            reference = "-"
            if alt != "-":
                alt_new = alt[len(r) :]
                alternative = alt_new
            else:
                alternative = str(alt)
        else:
            return p, r, alt
    return position, reference, alternative


def check_min_req_GSvar(row):
    """
    checking the presence of mandatory columns
    :param row: dictionary of a GSvar row
    :return: boolean, True if min req met
    """
    if (
        "#chr" in row.keys()
        and "start" in row.keys()
        and "end" in row.keys()
        and "ref" in row.keys()
        and "obs" in row.keys()
        and (
            "coding_and_splicing_details" in row.keys() or "coding" in row.keys() or "coding_and_splicing" in row.keys()
        )
    ):
        return True
    return False


def determine_variant_type(record, alternative):
    vt = VariationType.UNKNOWN
    if record.is_snp:
        vt = VariationType.SNP
    elif record.is_indel:
        if abs(len(alternative) - len(record.REF)) % 3 == 0:  # no frameshift
            if record.is_deletion:
                vt = VariationType.DEL
            else:
                vt = VariationType.INS
        else:  # frameshift
            if record.is_deletion:
                vt = VariationType.FSDEL
            else:
                vt = VariationType.FSINS
    return vt


def determine_zygosity(record):
    genotye_dict = {"het": False, "hom": True, "ref": True}
    isHomozygous = False
    if "HOM" in record.INFO:
        isHomozygous = record.INFO["HOM"] == 1
    elif "SGT" in record.INFO:
        zygosity = record.INFO["SGT"].split("->")[1]
        if zygosity in genotye_dict:
            isHomozygous = genotye_dict[zygosity]
        else:
            if zygosity[0] == zygosity[1]:
                isHomozygous = True
            else:
                isHomozygous = False
    else:
        for sample in record.samples:
            if "GT" in sample.data:
                isHomozygous = sample.data["GT"] == "1/1"
    return isHomozygous


def read_GSvar(filename, pass_only=True):
    """
    reads GSvar and tsv files (tab sep files in context of genetic variants), omitting and warning about rows missing
    mandatory columns
    :param filename: /path/to/file
    :return: list epytope variants
    """
    global ID_SYSTEM_USED
    RE = re.compile("(\w+):([\w.]+):([&\w]+):\w*:exon(\d+)\D*\d*:(c.\D*([_\d]+)\D*):(p.\D*(\d+)\w*)")

    # list of mandatory (meta)data
    exclusion_list = [
        "start",
        "end",
        "#chr",
        "ref",
        "obs",
        "gene",
        "tumour_genotype",
        "coding_and_splicing_details",
        "variant_details",
        "variant_type",
        "coding_and_splicing",
    ]

    list_vars = list()
    lines = list()
    transcript_ids = []
    dict_vars = {}

    cases = 0

    with open(filename, "rt") as tsvfile:
        tsvreader = csv.DictReader((row for row in tsvfile if not row.startswith("##")), delimiter="\t")
        for row in tsvreader:
            if not check_min_req_GSvar(row):
                logger.warning("read_GSvar: Omitted row! Mandatory columns not present in: \n" + str(row) + ".")
                continue
            lines.append(row)

    # get list of additional metadata
    metadata_list = set(tsvreader.fieldnames) - set(exclusion_list)

    for mut_id, line in enumerate(lines):
        if "filter" in line and pass_only and line["filter"].strip():
            continue
        genome_start = int(line["start"]) - 1
        genome_stop = int(line["end"]) - 1
        chrom = line["#chr"]
        ref = line["ref"]
        alt = line["obs"]
        gene = line.get("gene", "")

        isHomozygous = (
            True
            if (
                ("tumour_genotype" in line)
                and (line["tumour_genotype"].split("/")[0] == line["tumour_genotype"].split("/")[1])
            )
            else False
        )

        # old GSvar version
        if "coding_and_splicing_details" in line:
            mut_type = line.get("variant_details", "")
            annots = RE.findall(line["coding_and_splicing_details"])
        else:
            mut_type = line.get("variant_type", "")
            # Gene, transcript number, type, impact, exon/intron number, HGVS.c, HGVS.p, Pfam
            annots = RE.findall(line["coding_and_splicing"])
        isyn = mut_type == "synonymous_variant"

        """
        Enum for variation types:
        type.SNP, type.DEL, type.INS, type.FSDEL, type.FSINS, type.UNKNOWN
        """
        vt = VariationType.UNKNOWN
        if mut_type == "missense_variant" or "missense_variant" in mut_type:
            vt = VariationType.SNP
        elif mut_type == "frameshift_variant":
            if (ref == "-") or (len(ref) < len(alt)):
                vt = VariationType.FSINS
            else:
                vt = VariationType.FSDEL
        elif mut_type == "inframe_deletion":
            vt = VariationType.DEL
        elif mut_type == "inframe_insertion":
            vt = VariationType.INS

        coding = dict()

        for annot in annots:
            a_gene, transcript_id, a_mut_type, exon, trans_coding, trans_pos, prot_coding, prot_start = annot
            if "NM" in transcript_id:
                ID_SYSTEM_USED = EIdentifierTypes.REFSEQ
            if "stop_gained" not in mut_type:
                if not gene:
                    gene = a_gene
                if not mut_type:
                    mut_type = a_mut_type

                # TODO with the next epytope release we can deal with transcript id version
                transcript_id = transcript_id.split(".")[0]

                coding[transcript_id] = MutationSyntax(
                    transcript_id, int(trans_pos.split("_")[0]) - 1, int(prot_start) - 1, trans_coding, prot_coding
                )
                transcript_ids.append(transcript_id)
        if coding:
            var = Variant(
                mut_id,
                vt,
                chrom.strip("chr"),
                int(genome_start),
                ref.upper(),
                alt.upper(),
                coding,
                isHomozygous,
                isSynonymous=isyn,
            )
            var.gene = gene

            # metadata logging
            for meta_name in metadata_list:
                var.log_metadata(meta_name, line.get(meta_name, ""))

            dict_vars[var] = var
            list_vars.append(var)

    transToVar = {}

    # fix because of memory/timing issues due to combinatorial explosion
    for variant in list_vars:
        for trans_id in variant.coding.keys():
            transToVar.setdefault(trans_id, []).append(variant)

    for tId, vs in transToVar.items():
        if len(vs) > 10:
            cases += 1
            for v in vs:
                vs_new = Variant(v.id, v.type, v.chrom, v.genomePos, v.ref, v.obs, v.coding, True, v.isSynonymous)
                vs_new.gene = v.gene
                for m in metadata_list:
                    vs_new.log_metadata(m, v.get_metadata(m)[0])
                dict_vars[v] = vs_new
    return dict_vars.values(), transcript_ids, metadata_list


def read_vcf(filename, pass_only=True):
    """
    reads vcf files
    returns a list of epytope variants
    :param filename: /path/to/file
    :param boolean pass_only: only consider variants that passed the filter (default: True)
    :return: list of epytope variants
    """
    global ID_SYSTEM_USED

    SNPEFF_KEY = "ANN"
    MISSENSE_CLASS = "missense_variant"

    variants = list()
    with open(filename, "rt") as tsvfile:
        vcf_reader = vcf.Reader(tsvfile)
        variants = [r for r in vcf_reader]

    # get lists of additional metadata
    metadata_list = set(vcf_reader.infos.keys())
    format_list = set(vcf_reader.formats.keys())
    final_metadata_list = []

    dict_vars = {}
    list_vars = []
    transcript_ids = []

    for num, record in enumerate(variants):
        chromosome = record.CHROM.strip("chr")
        genomic_position = record.POS
        variation_dbid = record.ID
        reference = str(record.REF)
        alternative_list = record.ALT
        record_filter = record.FILTER

        if pass_only and record_filter:
            continue

        """
        Enum for variation types:
        type.SNP, type.DEL, type.INS, type.FSDEL, type.FSINS, type.UNKNOWN

        VARIANT INCORP IN EPYTOPE

        SNP => seq[pos] = OBS (replace)
        INSERTION => seqp[pos:pos] = obs (insert at that position)
        DELETION => s = slice(pos, pos+len(ref)) (create slice that will be removed)
		            del seq[s] (remove)
        """
        for alt in list(filter(None, alternative_list)):
            isHomozygous = determine_zygosity(record)
            vt = determine_variant_type(record, alt)
            if record.INFO.get(SNPEFF_KEY, False):
                isSynonymous = False
                coding = dict()
                types = []
                for annraw in record.INFO[SNPEFF_KEY]:
                    annots = annraw.split("|")
                    if len(annots) != 16:
                        logger.warning(
                            "read_vcf: Omitted row! Mandatory columns not present in annotation field (ANN). \n Have you annotated your VCF file with SnpEff?"
                        )
                        continue
                    (
                        obs,
                        a_mut_type,
                        impact,
                        a_gene,
                        a_gene_id,
                        feature_type,
                        transcript_id,
                        exon,
                        tot_exon,
                        trans_coding,
                        prot_coding,
                        cdna,
                        cds,
                        aa,
                        distance,
                        warnings
                    ) = annots

                    if MISSENSE_CLASS in annots: 
                        
                        types.append(a_mut_type)
                        tpos = 0
                        ppos = 0
                        positions = ""
                        isSynonymous = a_mut_type == "synonymous_variant"
                        gene = a_gene

                        # get cds/protein positions and convert mutation syntax to epytope format
                        if trans_coding != "":
                            positions = re.findall(r"\d+", trans_coding)
                            ppos = int(positions[0]) - 1

                        if prot_coding != "":
                            positions = re.findall(r"\d+", prot_coding)
                            tpos = int(positions[0]) - 1

                        # TODO with the new epytope release we will support transcript IDs with version
                        transcript_id = transcript_id.split(".")[0]
                        if "NM" in transcript_id:
                            ID_SYSTEM_USED = EIdentifierTypes.REFSEQ
                        
                        # take online coding variants into account, epytope cannot deal with stop gain variants right now
                        if not prot_coding or "stop_gained" in a_mut_type:
                            continue
                        
                        coding[transcript_id] = MutationSyntax(transcript_id, ppos, tpos, trans_coding, prot_coding)
                        transcript_ids.append(transcript_id)

            else:
                print('No ANN annotation. Skipping...')
                continue

            if coding:
                pos, reference, alternative = get_epytope_annotation(vt, genomic_position, reference, str(alt))
                var = Variant(
                    "line" + str(num),
                    vt,
                    chromosome,
                    pos,
                    reference,
                    alternative,
                    coding,
                    isHomozygous,
                    isSynonymous,
                )
                var.gene = gene
                var.log_metadata("vardbid", variation_dbid)
                final_metadata_list.append("vardbid")
                for metadata_name in metadata_list:
                    if metadata_name in record.INFO:
                        final_metadata_list.append(metadata_name)
                        var.log_metadata(metadata_name, record.INFO[metadata_name])
                for sample in record.samples:
                    for format_key in format_list:
                        if getattr(sample.data, format_key, None) is None:
                            logger.warning(
                                "FORMAT entry {entry} not defined for {genotype}. Skipping.".format(
                                    entry=format_key, genotype=sample.sample
                                )
                            )
                            continue
                        format_header = "{}.{}".format(sample.sample, format_key)
                        final_metadata_list.append(format_header)
                        if isinstance(sample[format_key], list):
                            format_value = ",".join([str(i) for i in sample[format_key]])
                        else:
                            format_value = sample[format_key]
                        var.log_metadata(format_header, format_value)
                dict_vars[var] = var
                list_vars.append(var)
            else:
                logger.error("No supported variant annotation string found. Skipping...")
                continue

            transToVar = {}
    
    # fix because of memory/timing issues due to combinatorial explosion
    for variant in list_vars:
        for trans_id in variant.coding.keys():
            transToVar.setdefault(trans_id, []).append(variant)
            
    for tId, vs in transToVar.items():
        if len(vs) > 10:
            for v in vs:
                vs_new = Variant(v.id, v.type, v.chrom, v.genomePos, v.ref, v.obs, v.coding, True, v.isSynonymous)
                vs_new.gene = v.gene
                for m in metadata_name:
                    vs_new.log_metadata(m, v.get_metadata(m))
                dict_vars[v] = vs_new

    return dict_vars.values(), transcript_ids, final_metadata_list


def read_peptide_input(filename):
    peptides = []
    metadata = []

    """expected columns (min required): id sequence"""
    with open(filename, "r") as peptide_input:
        # enable listing of protein names for each peptide
        csv.field_size_limit(600000)
        reader = csv.DictReader(peptide_input, delimiter="\t")
        for row in reader:
            pep = Peptide(row["sequence"])

            for col in row:
                if col != "sequence":
                    pep.log_metadata(col, row[col])
                    metadata.append(col)
            peptides.append(pep)

    metadata = set(metadata)
    return peptides, metadata


# parse protein_groups of MaxQuant output to get protein intensity values
def read_protein_quant(filename):
    # protein id: sample1: intensity, sample2: intensity:
    intensities = {}

    with open(filename, "r") as inp:
        inpreader = csv.DictReader(inp, delimiter="\t")
        for row in inpreader:
            if "REV" in row["Protein IDs"]:
                pass
            else:
                valuedict = {}
                for key, val in row.iteritems():
                    if "LFQ intensity" in key:
                        valuedict[key.replace("LFQ intensity ", "").split("/")[-1]] = val
                for p in row["Protein IDs"].split(";"):
                    if "sp" in p:
                        intensities[p.split("|")[1]] = valuedict
    return intensities


# parse different expression analysis results (DESeq2), link log2fold changes to transcripts/genes
def read_diff_expression_values(filename):
    # feature id: log2fold changes
    fold_changes = {}

    with open(filename, "r") as inp:
        inp.readline()
        for row in inp:
            values = row.strip().split("\t")
            fold_changes[values[0]] = values[1]

    return fold_changes


# parse ligandomics ID output, peptide sequences, scores and median intensity
def read_lig_ID_values(filename):
    # sequence: score median intensity
    intensities = {}

    with open(filename, "r") as inp:
        reader = csv.DictReader(inp, delimiter=",")
        for row in reader:
            seq = re.sub("[\(].*?[\)]", "", row["sequence"])
            intensities[seq] = (row["fdr"], row["intensity"])

    return intensities


def create_protein_column_value(pep):

    all_proteins = [
        transcriptProteinMap[transcript.transcript_id.split(":")[0]] for transcript in set(pep.get_all_transcripts())
    ]
    return ",".join(set([item for sublist in all_proteins for item in sublist]))


def create_transcript_column_value(pep):
    return ",".join(set([transcript.transcript_id.split(":")[0] for transcript in set(pep.get_all_transcripts())]))


def create_mutationsyntax_column_value(pep, pep_dictionary):
    syntaxes = []
    for variant in set(pep_dictionary[pep]):
        for coding in variant.coding:
            syntaxes.append(variant.coding[coding])
    return ",".join(set([mutationSyntax.aaMutationSyntax for mutationSyntax in syntaxes]))


def create_mutationsyntax_genome_column_value(pep, pep_dictionary):
    syntaxes = []
    for variant in set(pep_dictionary[pep]):
        for coding in variant.coding:
            syntaxes.append(variant.coding[coding])
    return ",".join(set([mutationSyntax.cdsMutationSyntax for mutationSyntax in syntaxes]))


def create_gene_column_value(pep, pep_dictionary):
    return ",".join(set([variant.gene for variant in set(pep_dictionary[pep])]))

def create_reference_column_value(pep, ref):
    return True if sum(aa1 != aa2 for aa1, aa2 in zip(pep, ref)) <= 1 and len(pep) == len(ref) else False
    

def create_variant_pos_column_value(pep, pep_dictionary):
    return ",".join(set(["{}".format(variant.genomePos) for variant in set(pep_dictionary[pep])]))


def create_variant_chr_column_value(pep, pep_dictionary):
    return ",".join(set(["{}".format(variant.chrom) for variant in set(pep_dictionary[pep])]))


def create_variant_type_column_value(pep, pep_dictionary):
    types = {0: "SNP", 1: "DEL", 2: "INS", 3: "FSDEL", 4: "FSINS", 5: "UNKNOWN"}
    return ",".join(set([types[variant.type] for variant in set(pep_dictionary[pep])]))


def create_variant_syn_column_value(pep, pep_dictionary):
    return ",".join(set([str(variant.isSynonymous) for variant in set(pep_dictionary[pep])]))


def create_variant_hom_column_value(pep, pep_dictionary):
    return ",".join(set([str(variant.isHomozygous) for variant in set(pep_dictionary[pep])]))


def create_coding_column_value(pep, pep_dictionary):
    return ",".join(set([str(variant.coding) for variant in set(pep_dictionary[pep])]))


def create_metadata_column_value(pep, c, pep_dictionary):
    meta = set(
        [
            str(variant.get_metadata(c)[0])
            for variant in set(pep_dictionary[pep[0]])
            if len(variant.get_metadata(c)) != 0
        ]
    )
    if len(meta) is 0:
        return np.nan
    else:
        return ",".join(meta)


def create_wt_seq_column_value(pep, wtseqs):
    transcripts = [transcript for transcript in set(pep["sequence"].get_all_transcripts())]
    wild_type = set(
        [
            str(wtseqs["{}_{}".format(str(pep["sequence"]), transcript.transcript_id)])
            for transcript in transcripts
            if bool(transcript.vars) and "{}_{}".format(str(pep["sequence"]), transcript.transcript_id) in wtseqs
        ]
    )
    if len(wild_type) is 0:
        return np.nan
    else:
        return ",".join(wild_type)


def create_quant_column_value(row, dict):
    if row[1] in dict:
        value = dict[row[1]]
    else:
        value = np.nan
    return value


# defined as : RPKM = (10^9 * C)/(N * L)
# L = exon length in base-pairs for a gene
# C = Number of reads mapped to a gene in a single sample
# N = total (unique)mapped reads in the sample
def create_expression_column_value_for_result(row, dict, deseq, gene_id_lengths):
    ts = row["gene"].split(",")
    values = []
    if deseq:
        for t in ts:
            if t in dict:
                values.append(dict[t])
            else:
                values.append(np.nan)
    else:
        for t in ts:
            if t in dict:
                if t in gene_id_lengths:
                    values.append(
                        (10.0**9 * float(dict[t]))
                        / (
                            float(gene_id_lengths[t])
                            * sum(
                                [
                                    float(dict[k])
                                    for k in dict.keys()
                                    if ((not k.startswith("__")) & (k in gene_id_lengths))
                                ]
                            )
                        )
                    )
                else:
                    values.append(
                        (10.0**9 * float(dict[t]))
                        / (
                            float(len(row[0].get_all_transcripts()[0]))
                            * sum(
                                [
                                    float(dict[k])
                                    for k in dict.keys()
                                    if ((not k.startswith("__")) & (k in gene_id_lengths))
                                ]
                            )
                        )
                    )
                    logger.warning(
                        "FKPM value will be based on transcript length for {gene}. Because gene could not be found in the DB".format(
                            gene=t
                        )
                    )
            else:
                values.append(np.nan)
    values = ["{0:.2f}".format(value) for value in values]
    return ",".join(values)


def create_quant_column_value_for_result(row, dict, swissProtDict, key):
    all_proteins = [swissProtDict[x.transcript_id.split(":")[0]] for x in set(row[0].get_all_transcripts())]
    all_proteins_filtered = set([item for sublist in all_proteins for item in sublist])
    values = []
    for p in all_proteins_filtered:
        if p in dict:
            if int(dict[p][key]) > 0:
                values.append(math.log(int(dict[p][key]), 2))
            else:
                values.append(int(dict[p][key]))
    if len(values) is 0:
        return np.nan
    else:
        return ",".join(set([str(v) for v in values]))


def create_ligandomics_column_value_for_result(row, lig_id, val, wild_type):
    if wild_type:
        seq = row["wt sequence"]
    else:
        seq = row["sequence"]
    if seq in lig_id:
        return lig_id[seq][val]
    else:
        return ""


def get_protein_ids_for_transcripts(idtype, transcripts, ensembl_url, reference):
    result = {}
    result_swissProt = {}

    biomart_url = "{}/biomart/martservice?query=".format(ensembl_url)
    
    biomart_head = """
    <?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE Query>
        <Query client="true" processor="TSV" limit="-1" header="1" uniqueRows = "1" >
            <Dataset name="%s" config="%s">
    """.strip()
    biomart_tail = """
            </Dataset>
        </Query>
    """.strip()
    biomart_filter = """<Filter name="%s" value="%s" filter_list=""/>"""
    biomart_attribute = """<Attribute name="%s"/>"""

    ENSEMBL = False
    if idtype == EIdentifierTypes.ENSEMBL:
        idname = "ensembl_transcript_id"
        ENSEMBL = True
    elif idtype == EIdentifierTypes.REFSEQ:
        idname = "refseq_mrna"

    input_lists = []

    # too long requests will fail
    if len(transcripts) > 200:
        input_lists = [transcripts[i : i + 3] for i in range(0, len(transcripts), 3)]

    else:
        input_lists += [transcripts]

    attribut_swissprot = "uniprot_swissprot_accession" if reference == "GRCh37" else "uniprotswissprot"

    tsvselect = []
    for l in input_lists:
        rq_n = (
            biomart_head % ("hsapiens_gene_ensembl", "default")
            + biomart_filter % (idname, ",".join(l))
            + biomart_attribute % ("ensembl_peptide_id")
            + biomart_attribute % (attribut_swissprot)
            + biomart_attribute % ("refseq_peptide")
            + biomart_attribute % (idname)
            + biomart_tail
        )

        # DictReader returns byte object that is transformed into a string by '.decode('utf-8')'
        tsvreader = csv.DictReader(
            urllib.request.urlopen(biomart_url + urllib.parse.quote(rq_n)).read().decode("utf-8").splitlines(),
            dialect="excel-tab",
        )

        tsvselect += [x for x in tsvreader]

    swissProtKey = "UniProt/SwissProt Accession" if reference == "GRCh37" else "UniProtKB/Swiss-Prot ID"

    if ENSEMBL:
        key = "Ensembl Transcript ID" if reference == "GRCh37" else "Transcript stable ID"
        protein_key = "Ensembl Protein ID" if reference == "GRCh37" else "Protein stable ID"
        for dic in tsvselect:
            if dic[key] in result:
                merged = result[dic[key]] + [dic[protein_key]]
                merged_swissProt = result_swissProt[dic[key]] + [dic[swissProtKey]]
                result[dic[key]] = merged
                result_swissProt[dic[key]] = merged_swissProt
            else:
                result[dic[key]] = [dic[protein_key]]
                result_swissProt[dic[key]] = [dic[swissProtKey]]
    else:
        key = "RefSeq mRNA [e.g. NM_001195597]"
        for dic in tsvselect:
            if dic[key] in result:
                merged = result[dic[key]] + [dic["RefSeq Protein ID [e.g. NP_001005353]"]]
                merged_swissProt = result_swissProt[dic[key]] + [dic[swissProtKey]]
                result[dic[key]] = merged
                result_swissProt[dic[key]] = merged_swissProt
            else:
                result[dic[key]] = [dic["RefSeq Protein ID [e.g. NP_001005353]"]]
                result_swissProt[dic[key]] = [dic[swissProtKey]]

    return result, result_swissProt


def get_matrix_max_score(allele, length):
    allele_model = "%s_%i" % (allele, length)
    try:
        pssm = getattr(
            __import__("epytope.Data.pssms.syfpeithi.mat." + allele_model, fromlist=[allele_model]), allele_model
        )
        return sum([max(scrs.values()) for pos, scrs in pssm.items()])
    except:
        return np.nan


def create_affinity_values(allele, length, j, method, max_scores, allele_strings):
    if not pd.isnull(j):
        if "syf" in method:
            return max(
                0, round((100.0 / float(max_scores[allele_strings[("%s_%s" % (str(allele), length))]]) * float(j)), 2)
            )
        else:
            # convert given affinity score in range [0,1] back to IC50 affinity value
            return round((50000 ** (1.0 - float(j))), 2)
    else:
        return np.nan


def create_binder_values(pred_value, method, thresholds):
    if not pd.isnull(pred_value):
        if "syf" in method:
            return True if pred_value > thresholds[method] else False
        else:
            return True if pred_value <= thresholds[method.lower()] else False
    else:
        return np.nan


def generate_wt_seqs(peptides):
    wt_dict = {}

    r = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
    d_pattern = re.compile("([a-zA-Z]+)([0-9]+)")
    for x in peptides:
        trans = x.get_all_transcripts()
        for t in trans:
            mut_seq = [a for a in x]
            protein_pos = x.get_protein_positions(t.transcript_id)
            not_available = False
            variant_available = False
            for p in protein_pos:
                variant_dic = x.get_variants_by_protein_position(t.transcript_id, p)
                variant_available = bool(variant_dic)
                for key in variant_dic:
                    var_list = variant_dic[key]
                    for v in var_list:
                        mut_syntax = v.coding[t.transcript_id.split(":")[0]].aaMutationSyntax
                        if v.type in [3, 4, 5] or "?" in mut_syntax:
                            not_available = True
                        elif v.type in [1]:
                            m = d_pattern.match(mut_syntax.split(".")[1])
                            wt = SeqUtils.seq1(m.groups()[0])
                            mut_seq.insert(key, wt)
                        elif v.type in [2]:
                            not_available = True
                        else:
                            m = r.match(mut_syntax.split(".")[1])
                            if m is None:
                                not_available = True
                            else:
                                wt = SeqUtils.seq1(m.groups()[0])
                                mut_seq[key] = wt
            if not_available:
                wt_dict["{}_{}".format(str(x), t.transcript_id)] = np.nan
            elif variant_available:
                wt_dict["{}_{}".format(str(x), t.transcript_id)] = "".join(mut_seq)
    return wt_dict


# TODO potential improvement in epytope
def create_peptide_variant_dictionary(peptides):
    pep_to_variants = {}
    for pep in peptides:
        transcript_ids = [x.transcript_id for x in set(pep.get_all_transcripts())]
        variants = []
        for t in transcript_ids:
            variants.extend([v for v in pep.get_variants_by_protein(t)])
        pep_to_variants[pep] = variants
    return pep_to_variants


# TODO replace by epytope function once released
def is_created_by_variant(peptide):
    transcript_ids = [x.transcript_id for x in set(peptide.get_all_transcripts())]
    for t in transcript_ids:
        p = peptide.proteins[t]
        varmap = p.vars
        for pos, vars in varmap.items():
            for var in vars:
                if var.type in [VariationType.FSDEL, VariationType.FSINS]:
                    if peptide.proteinPos[t][0] + len(peptide) > pos:
                        return True
                else:
                    for start_pos in peptide.proteinPos[t]:
                        positions = list(range(start_pos, start_pos + len(peptide)))
                        if pos in positions:
                            return True
    return False

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO
import pandas as pd
import logging
import warnings


from epytope.Core.Base import MetadataLogger
from epytope.Core.Transcript import Transcript
import epytope.Core.Generator as generator
from epytope.Core.Variant import Variant, VariationType, MutationSyntax
from epytope.Core.Generator import _incorp

EAdapterFields = (lambda **enums: type('Enum', (), enums))(GENE=0, STRAND=1, SEQ=2, TRANSID=3, PROTID=4)
EIdentifierTypes = (lambda **enums: type('Enum', (), enums))(ENSEMBL=0, REFSEQ=1, PREDREFSEQ=2, UNIPROT=3, GENENAME=4, HGNC=4)
#VariationType = (lambda **enums: type('Enum', (), enums))(SNP=0, DEL=1, INS=2, FSDEL=3, FSINS=4, UNKNOWN=5)


# _incorp = {
#             VariationType.DEL: generator._incorp_deletion,
#             VariationType.FSDEL: generator._incorp_deletion,
#             VariationType.INS: generator._incorp_insertion,
#             VariationType.FSINS: generator._incorp_insertion,
#             VariationType.SNP: generator._incorp_snp}
         

REVERS = -1

class GenbankRetrievalError(Exception):
    pass

def fetch_genbank_record(accession_number):
    """
    Retrieve GenBank record for a given accession number from NCBI.

    Parameters:
    - accession_number (str): The accession number of the sequence.

    Returns:
    - str: The GenBank record in text format.

    Raises:
    - ValueError: If the accession number is not provided.
    - GenbankRetrievalError: If the retrieval fails (e.g., invalid accession number or network issues).
    """

    # Validate input
    if not accession_number:
        raise ValueError("Accession number must be provided.")

    try:
        # Set your email address (required for accessing NCBI data)
        Entrez.email = "your.email@example.com"

        # Fetch the GenBank record for the specified accession number
        genbank = Entrez.efetch(db="nuccore", id=accession_number, rettype="gb", retmode="text")
        handle= SeqIO.read(genbank, "genbank")
        genbank.close()
    
        return handle

    except Entrez.HTTPError as e:
        # Handle HTTP errors (e.g., 404 Not Found, 500 Internal Server Error)
        raise GenbankRetrievalError(f"HTTP Error: {str(e)}")
    except Exception as e:
        # Handle other exceptions
        raise GenbankRetrievalError(f"Failed to retrieve GenBank record: {str(e)}")

def fetch_transcripts(locus_tags, seq_record, intergene_length=1):
    cds_list_plus = []
    cds_list_minus = []
    
    try:
        # Replace "your.email@example.com" with your actual email address
                
        intergenic_records=[]
        #print(StringIO(transcripts))
        for feature in seq_record.features:
            tag = feature.qualifiers.get('locus_tag')
            geneid = feature.qualifiers.get('gene')
            proteinid = feature.qualifiers.get('protein_id')
            db_xref = feature.qualifiers.get('db_xref')

            if tag:
                if tag[0] in locus_tags:
                    if feature.type == 'CDS' and "translation" in feature.qualifiers:                              
                        translation = feature.qualifiers["translation"][0]
                        nuc_seq = feature.location.extract(seq_record).seq
                        prot_seq = nuc_seq.translate()

                        if prot_seq.replace('*', '') == translation:
                            #mystart = feature.location._start.position
                            #myend = feature.location._end.position

                            if feature.strand == -1:
                                cds_list_minus.append(feature.location)

                            elif feature.strand == 1:
                                cds_list_plus.append((feature.location))

                            if len(nuc_seq) >= intergene_length:
                                
                                intergenic_records =    {EAdapterFields.SEQ: nuc_seq,
                                                         EAdapterFields.GENE: geneid,
                                                         EAdapterFields.STRAND: feature.strand,
                                                         EAdapterFields.TRANSID: tag[0],
                                                         EAdapterFields.PROTID: proteinid}
                                return intergenic_records
        return 'Sem intergene'
        #return intergenic_records
        #print("Transcripts:")
        #for transcript in transcripts:
        #    print(transcript)

    except GenbankRetrievalError as gre:
        print(f"Genbank Retrieval Error: {str(gre)}")

def generate_transcripts_from_variants_gb(variants):
    """
    Generates all possible transcript :class:`~epytope.Core.Transcript.Transcript` based on the given
    :class:`~epytope.Core.Variant.Variant`.

    The result is a generator.

    :param vars: A list of variants for which transcripts should be build
    :type vars: list(:class:`~epytope.Core.Variant.Variant`)
    :param: dbadapter: a DBAdapter to fetch the transcript sequences
    :type dbadapter: class:`~epytope.IO.ADBAdapter.ADBAdapter`
    :param id_type: The type of the transcript IDs used in annotation of variants (e.g. REFSEQ, ENSAMBLE)
    :type id_type: :func:`~epytope.IO.ADBAdapter.EIdentifierTypes`
    :return: A generator of transcripts with all possible variations determined by the given variant list
    :rtype: Generator(:class:`~epytope.Core.Transcript.Transcript)
    :invariant: Variants are considered to be annotated from forward strand, regardless of the transcripts real
                orientation
    """
    
    def _generate_combinations(tId, vs, seq, usedVs, offset, isReverse=False):
        """
        recursive variant combination generator
        """
        transOff = generate_transcripts_from_variants_gb.transOff
        #print "TransOffset ", transOff, tId,usedVs
        if vs:
            v = vs.pop()
            if v.isHomozygous:
                pos = v.coding[tId].tranPos + offset
                usedVs[pos] = v
                offset = _incorp.get(v.type, lambda a, b, c, d, e, f: e)(seq, v, tId, pos, offset, isReverse)

                for s in _generate_combinations(tId, vs, seq, usedVs, offset, isReverse):
                    yield s
            else:
                vs_tmp = vs[:]
                tmp_seq = seq[:]
                tmp_usedVs = usedVs.copy()

                for s in _generate_combinations(tId, vs_tmp, tmp_seq, tmp_usedVs, offset, isReverse):
                    yield s

                # update the transcript variant id
                generate_transcripts_from_variants_gb.transOff += 1
                pos = v.coding[tId].tranPos + offset
                usedVs[pos] = v
                offset = _incorp.get(v.type, lambda a, b, c, d, e, f: e)(seq, v, tId, pos, offset, isReverse)
                for s in _generate_combinations(tId, vs, seq, usedVs, offset, isReverse):
                    yield s
        else:
            yield tId+":epytope_%i"%transOff, seq, usedVs

    #1) get all transcripts and sort the variants to transcripts

    #For a transcript do:
        #A) get transcript sequences
        #B) generate all possible combinations of variants
        #C) apply variants to transcript and generate transcript object

    transToVar = {}
    for v in variants:
       for trans_id in v.coding.keys():
            #TODO: verify if GU affect something
            if 'GU' in trans_id:
                transToVar.setdefault(trans_id, []).append(v)
    
    #for tId in variants:
    for tId, vs in transToVar.items():
        
        print('Variants Data')
        
        for i in vs:
            print(vars(i))

        Entrez.email = "your.email@example.com"
        seq_record = fetch_genbank_record('NC_045512.2')
        query = fetch_transcripts(tId, seq_record=seq_record)

        if query is None:
            logging.warning("Transcript with ID %s not found in DB or sequence unavailable"%tId)
            continue

        tSeq = query[EAdapterFields.SEQ]
        geneid = query[EAdapterFields.TRANSID]
        strand = query[EAdapterFields.STRAND]
        
        vs.sort(key=lambda v: v.genomePos-1
                if v.type in [VariationType.FSINS, VariationType.INS]
                else v.genomePos, reverse=True)
        if not generator._check_for_problematic_variants(vs):
            warnings.warn("Intersecting variants found for Transcript %s"%tId)
            continue
        generate_transcripts_from_variants_gb.transOff = 0
        varRef = None
        for tId, varSeq, varComb in _generate_combinations(tId, vs, list(tSeq), {}, 0, isReverse=strand == REVERS):
            print('GENERATE COMBINATIONS  -------------------------------')
            print(tId)
            print("".join(varSeq))
            print(varComb)
            transcript_class = Transcript("".join(varSeq), geneid, tId, vars=varComb)
            transcript_class.reference = tSeq
            print(vars(transcript_class))
            print(' FINISH COMBINATIONS  -------------------------------')
            yield transcript_class

def generate_proteins_from_transcripts_gb(transcripts, table='Standard', stop_symbol='*', to_stop=True, cds=False):
        """
        Enables the translation from a :class:`~epytope.Core.Transcript.Transcript` to a
        :class:`~epytope.Core.Protein.Protein` instance. The result is a generator.

        The result is a generator.

        :param transcripts:  A list of or a single transcripts to translate
        :type transcripts: list(:class:`~epytope.Core.Transcript.Transcript`) or :class:`~epytope.Core.Transcript.Transcript`
        :param str table: Which codon table to use? This can be either a name (string), an NCBI identifier (integer),
                          or a CodonTable object (useful for non-standard genetic codes). Defaults to the 'Standard'
                          table
        :param str stop_symbol: Single character string, what to use for any terminators, defaults to the asterisk, '*'
        :param bool to_stop: Translates sequence and passes any stop codons if False (default True)(translated as the
                             specified stop_symbol). If True, translation is terminated at the first in frame stop
                             codon (and the stop_symbol is not appended to the returned protein sequence)
        :param bool cds: Boolean, indicates this is a complete CDS. If True, this checks the sequence starts with a
                         valid alternative start codon (which will be translated as methionine, M), that the sequence
                         length is a multiple of three, and that there is a single in frame stop codon at the end
                         (this will be excluded from the protein sequence, regardless of the to_stop option).
                         If these tests fail, an exception is raised
        :returns: The protein that corresponds to the transcript
        :rtype: Generator(:class:`~epytope.Core.Protein.Protein`)
        :raises ValueError: If incorrect table argument is pasted
        :raises TranslationError: If sequence is not multiple of three, or first codon is not a start codon, or last
                                  codon ist not a stop codon, or an extra stop codon was found in frame, or codon is
                                  non-valid

        """

        if isinstance(transcripts, Transcript):
            transcripts = [transcripts]

        for t in transcripts:
            if not isinstance(t, Transcript):
                raise ValueError("An element of specified input is not of type Transcript")
            # translate to a protein sequence
            #if len(str(self)) % 3 != 0:
            #    raise ValueError('ERROR while translating: lenght of transcript %s is no multiple of 3, the transcript is:\n %s' % (self.transcript_id, self))

            #TODO warn if intrasequence stops - biopython warns if  % 3 != 0
            prot_seq = str(t.translate(table=table, stop_symbol=stop_symbol, to_stop=to_stop, cds=cds))
            
            #Get reference
            refSeq = t.reference
            ref_prot_seq = str(refSeq.translate(table=table, stop_symbol=stop_symbol, to_stop=to_stop, cds=cds))

            new_vars = dict()
            for pos, var in t.vars.items():
                if not var.isSynonymous:
                    prot_pos = pos // 3
                    new_vars.setdefault(prot_pos, []).append(var)

            gene_id = t.gene_id

            protein_class = Protein(prot_seq, gene_id, t.transcript_id, t, new_vars)
            protein_class.reference = ref_prot_seq
            yield protein_class

def generate_peptides_from_proteins_gb(proteins, window_size, peptides=None):
    """
    Creates all :class:`~epytope.Core.Peptide.Peptide` for a given window size, from a given
    :class:`~epytope.Core.Protein.Protein`.

    The result is a generator.

    :param proteins: (Iterable of) protein(s) from which a list of unique peptides should be generated
    :type proteins: list(:class:`~epytope.Core.Protein.Protein`) or :class:`~epytope.Core.Protein.Protein`
    :param int window_size: Size of peptide fragments
    :param peptides: A list of peptides to update during peptide generation (usa case: Adding and updating Peptides of
                     newly generated Proteins)
    :type peptides: list(:class:`~epytope.Core.Peptide.Peptide`)
    :return: A unique generator of peptides
    :rtype: Generator(:class:`~epytope.Core.Peptide.Peptide`)
    """

    def gen_peptide_info(protein):
        # Generate peptide sequences and returns the sequence
        # #and start position within the protein
        res = []

        seq = str(protein)
        refseq = str(protein.reference)

        for i in range(len(protein)+1-window_size):
            # generate peptide fragment
            end = i+window_size
            pep_seq = seq[i:end]
            ref_pep_seq = refseq[i:end]

            res.append((pep_seq, ref_pep_seq, i))
        
        return res

    if isinstance(peptides, Peptide):
        peptides = [peptides]

    final_peptides = {}

    if peptides:
        for p in peptides:
            if not isinstance(p, Peptide):
                raise ValueError("Specified list of Peptides contain non peptide objects")
            final_peptides[str(p)] = p

    if isinstance(proteins, Protein):
        proteins = [proteins]

    for prot in proteins:
        if not isinstance(prot, Protein):
            raise ValueError("Input does contain non protein objects.")
        # generate all peptide sequences per protein:
        for (seq, ref, pos) in gen_peptide_info(prot):
            if all(a in _allowed_aas for a in seq.upper()):
                t_id = prot.transcript_id
                if seq not in final_peptides:
                    pep_seq = Peptide(seq)
                    pep_seq.reference = Peptide(ref)
                    final_peptides[seq] = pep_seq

                final_peptides[seq].proteins[t_id] = prot
                final_peptides[seq].proteinPos[t_id].append(pos)

                #print(final_peptides.values())

    return iter(final_peptides.values())

def make_predictions_from_variants(
    variants_all,
    methods,
    tool_thresholds,
    use_affinity_thresholds,
    alleles,
    minlength,
    maxlength,
    martsadapter,
    protein_db,
    identifier,
    metadata,
    transcriptProteinMap,
):
    # list for all peptides and filtered peptides
    all_peptides = []
    all_peptides_filtered = []

    # dictionaries for syfpeithi matrices max values and allele mapping
    max_values_matrices = {}
    allele_string_map = {}

    # list to hold dataframes for all predictions
    pred_dataframes = []
    prots = [
        p
        for p in generate_proteins_from_transcripts_gb(
            generate_transcripts_from_variants_gb(variants_all)
        )
    ]

    for peplen in range(minlength, maxlength):
        peptide_gen = generate_peptides_from_proteins_gb(prots, peplen)
        peptides_var = [x for x in peptide_gen]
        peptides = [p for p in peptides_var if is_created_by_variant(p)]

        # filter out self peptides
        selfies = [str(p) for p in peptides if protein_db.exists(str(p))]
        filtered_peptides = [p for p in peptides if str(p) not in selfies]

        all_peptides = all_peptides + peptides
        all_peptides_filtered = all_peptides_filtered + filtered_peptides
     
        results = []
        results_ref = []
        references = []
        #print(filtered_peptides)
        if len(filtered_peptides) > 0:
            for method, version in methods.items():
                
                references.extend([ref.reference for ref in filtered_peptides])

                try:
                    predictor = EpitopePredictorFactory(method, version=version)
                    results.extend([predictor.predict(filtered_peptides, alleles=alleles)])
                    results_ref.extend([predictor.predict(references, alleles=alleles)])

                except:
                    logger.warning(
                        "Prediction for length {length} and allele {allele} not possible with {method} version {version}.".format(
                            length=peplen, allele=",".join([str(a) for a in alleles]), method=method, version=version
                        )
                    )
        
        # merge dataframes for multiple predictors
        if len(results) > 1:
            df = results[0].merge_results(results[1:])
            df_ref = results_ref[0].merge_results(results_ref[1:])
        elif len(results) == 1:
            df = results[0]
            df_ref = results_ref[0]
        else:

            continue
        
        df = pd.concat(results)
        df_ref = pd.concat(results_ref)
      
        # create method index and remove it from multi-column
        df = df.stack(level=1)
        df_ref = df_ref.stack(level=1)

        # merge remaining multi-column Allele and ScoreType
        df.columns = df.columns.map("{0[0]} {0[1]}".format)
        df_ref.columns = df_ref.columns.map("{0[0]} {0[1]}".format)

        # reset index to have indices as columns
        df.reset_index(inplace=True)
        df = df.rename(columns={"Method": "method", "Peptides": "sequence"})
        
        
        df_ref.reset_index(inplace=True)
        df_ref = df_ref.rename(columns={"Method": "method", "Peptides": "refseq"})
        #df_ref.columns = [str(col) + '_refseq' if col != 'refseq' else col for col in df_ref.columns]
        #df = df.add_suffix('_refseq')

        for a in alleles:
            conv_allele = "%s_%s%s" % (a.locus, a.supertype, a.subtype)
            allele_string_map["%s_%s" % (a, peplen)] = "%s_%i" % (conv_allele, peplen)
            max_values_matrices["%s_%i" % (conv_allele, peplen)] = get_matrix_max_score(conv_allele, peplen)

        pep_to_variants = create_peptide_variant_dictionary(df["sequence"].tolist())

        df["length"] = df["sequence"].map(len)
        df["chr"] = df["sequence"].map(lambda x: create_variant_chr_column_value(x, pep_to_variants))
        df["pos"] = df["sequence"].map(lambda x: create_variant_pos_column_value(x, pep_to_variants))
        df["gene"] = df["sequence"].map(lambda x: create_gene_column_value(x, pep_to_variants))
        df["transcripts"] = df["sequence"].map(create_transcript_column_value)
        #print(df)
        #df["proteins"] = df["sequence"].map(create_protein_column_value)
        df["variant type"] = df["sequence"].map(lambda x: create_variant_type_column_value(x, pep_to_variants))
        df["synonymous"] = df["sequence"].map(lambda x: create_variant_syn_column_value(x, pep_to_variants))
        df["homozygous"] = df["sequence"].map(lambda x: create_variant_hom_column_value(x, pep_to_variants))
        df["variant details (genomic)"] = df["sequence"].map(
            lambda x: create_mutationsyntax_genome_column_value(x, pep_to_variants)
        )
        df["variant details (protein)"] = df["sequence"].map(
            lambda x: create_mutationsyntax_column_value(x, pep_to_variants)
        )
        #TODO verificar affinity and score
        #TODO calcular para refseq
        df['refseq'] = references
        df['single_mutation'] = df.index.map(
            lambda i: create_reference_column_value(df.at[i, 'sequence'], references[i]))

        df = pd.merge(df, df_ref, on="refseq", how="left", suffixes=('', '_refseq'))

        df.to_csv("ref_result.tsv")

        for c in df.columns:
            if ("HLA-" in str(c) or "H-2-" in str(c)) and "Score" in str(c):
                idx = df.columns.get_loc(c)
                allele = c.rstrip(" Score")

                add_str = ""
                if "refseq" in str(c):
                    allele = c.rstrip(" Score_refseq")
                    add_str = "_refseq"

                df[c] = df[c].round(4)
                df.insert(
                    idx + 1,
                    "%s affinity%s" % (allele, add_str),
                    df.apply(
                        lambda x: create_affinity_values(
                            allele, int(x["length"]), float(x[c]), x["method"], max_values_matrices, allele_string_map
                        ),
                        axis=1,
                    ),
                )

                df.insert(
                    idx + 2,
                    "%s binder%s" % (allele, add_str),
                    df.apply(
                        lambda x: create_binder_values(float(x["%s Rank" % allele]), x["method"], tool_thresholds)
                        if "netmhc" in x["method"] and not use_affinity_thresholds
                        else create_binder_values(float(x["%s affinity" % allele]), x["method"], tool_thresholds),
                        axis=1,
                    ),
                )

                alleles_binders = df.columns[df.columns.str.contains(f'{allele} binder', regex = False)]
                if len(alleles_binders) == 2:
                    
                    df.insert(
                        idx + 3,
                        f"{allele} binder_changes" ,
                        df[f'{allele} binder'] != df[f'{allele} binder_refseq'])          


        df.columns = df.columns.str.replace("Score", "score")
        df.columns = df.columns.str.replace("Rank", "rank")

        for col in set(metadata):
            df[col] = df.apply(lambda row: create_metadata_column_value(row, col, pep_to_variants), axis=1)

        pred_dataframes.append(df)

    statistics = {
        "prediction_methods": [method + "-" + version for method, version in methods.items()],
        "number_of_variants": len(variants_all),
        "number_of_unique_peptides": [str(p) for p in all_peptides],
        "number_of_unique_peptides_after_filtering": [str(p) for p in all_peptides_filtered],
    }

    return pred_dataframes, statistics, all_peptides_filtered, prots


def make_predictions_from_peptides(
    peptides, methods, tool_thresholds, use_affinity_thresholds, alleles, protein_db, identifier, metadata
):
    # dictionaries for syfpeithi matrices max values and allele mapping
    max_values_matrices = {}
    allele_string_map = {}

    # list to hold dataframes for all predictions
    pred_dataframes = []

    # filter out self peptides if specified
    selfies = [str(p) for p in peptides if protein_db.exists(str(p))]
    peptides_filtered = [p for p in peptides if str(p) not in selfies]

    # sort peptides by length (for predictions)
    sorted_peptides = {}

    for p in peptides_filtered:
        length = len(str(p))
        if length in sorted_peptides:
            sorted_peptides[length].append(p)
        else:
            sorted_peptides[length] = [p]

    for peplen in sorted_peptides:
        all_peptides_filtered = sorted_peptides[peplen]
        results = []
        for method, version in methods.items():
            try:
                predictor = EpitopePredictorFactory(method, version=version)
                results.extend([predictor.predict(all_peptides_filtered, alleles=alleles)])
            except:
                logger.warning(
                    "Prediction for length {length} and allele {allele} not possible with {method} version {version}. No model available.".format(
                        length=peplen, allele=",".join([str(a) for a in alleles]), method=method, version=version
                    )
                )

        # merge dataframes for multiple predictors
        if len(results) > 1:
            df = results[0].merge_results(results[1:])
        elif len(results) == 1:
            df = results[0]
        else:
            continue

        # create method index and remove it from multi-column
        df = df.stack(level=1)

        # merge remaining multi-column Allele and ScoreType
        df.columns = df.columns.map("{0[0]} {0[1]}".format)

        # reset index to have indices as columns
        df.reset_index(inplace=True)
        df = df.rename(columns={"Method": "method", "Peptides": "sequence"})

        # create column containing the peptide lengths
        df.insert(2, "length", df["sequence"].map(len))

        for a in alleles:
            conv_allele = "%s_%s%s" % (a.locus, a.supertype, a.subtype)
            allele_string_map["%s_%s" % (a, peplen)] = "%s_%i" % (conv_allele, peplen)
            max_values_matrices["%s_%i" % (conv_allele, peplen)] = get_matrix_max_score(conv_allele, peplen)

        mandatory_columns = [
            "chr",
            "pos",
            "gene",
            "transcripts",
            "proteins",
            "variant type",
            "synonymous",
            "homozygous",
            "variant details (genomic)",
            "variant details (protein)",
        ]

        for header in mandatory_columns:
            if header not in metadata:
                df[header] = np.nan
            else:
                df[header] = df.apply(lambda row: row[0].get_metadata(header)[0], axis=1)

        for c in list(set(metadata) - set(mandatory_columns)):
            df[c] = df.apply(lambda row: row[0].get_metadata(c)[0], axis=1)

        for c in df.columns:
            if ("HLA-" in str(c) or "H-2-" in str(c)) and "Score" in str(c):
                idx = df.columns.get_loc(c)
                allele = c.rstrip(" Score")
                df[c] = df[c].round(4)
                df.insert(
                    idx + 1,
                    "%s affinity" % allele,
                    df.apply(
                        lambda x: create_affinity_values(
                            allele, int(x["length"]), float(x[c]), x["method"], max_values_matrices, allele_string_map
                        ),
                        axis=1,
                    ),
                )
                df.insert(
                    idx + 2,
                    "%s binder" % allele,
                    df.apply(
                        lambda x: create_binder_values(float(x["%s Rank" % allele]), x["method"], tool_thresholds)
                        if "netmhc" in x["method"] and not use_affinity_thresholds
                        else create_binder_values(float(x["%s affinity" % allele]), x["method"], tool_thresholds),
                        axis=1,
                    ),
                )

        df.columns = df.columns.str.replace("Score", "score")
        df.columns = df.columns.str.replace("Rank", "rank")

        pred_dataframes.append(df)

    # write prediction statistics
    statistics = {
        "prediction_methods": [method + "-" + version for method, version in methods.items()],
        "number_of_variants": 0,
        "number_of_unique_peptides": [str(p) for p in peptides],
        "number_of_unique_peptides_after_filtering": [str(p) for p in peptides_filtered],
    }
    return pred_dataframes, statistics


def __main__():
    parser = argparse.ArgumentParser(
        description="""EPAA - Epitope Prediction And Annotation \n Pipeline for prediction of MHC class I and II epitopes from variants or peptides for a list of specified alleles.
        Additionally predicted epitopes can be annotated with protein quantification values for the corresponding proteins, identified ligands, or differential expression values for the corresponding transcripts."""
    )
    parser.add_argument("-s", "--somatic_mutations", help="Somatic variants")
    parser.add_argument("-g", "--germline_mutations", help="Germline variants")
    parser.add_argument("-i", "--identifier", help="Dataset identifier")
    parser.add_argument("-p", "--peptides", help="File with one peptide per line")
    parser.add_argument("-l", "--max_length", help="Maximum peptide length")
    parser.add_argument("-ml", "--min_length", help="Minimum peptide length")
    parser.add_argument("-t", "--tools", help="Tools used for peptide predictions", required=True, type=str)
    parser.add_argument(
        "-tt",
        "--tool_thresholds",
        help="Customize thresholds of given tools using a json file",
        required=False,
        type=str,
    )
    parser.add_argument(
        "-at",
        "--use_affinity_thresholds",
        help="Use affinity instead of rank for thresholding",
        required=False,
        action="store_true",
    )
    parser.add_argument("-sv", "--versions", help="File containing parsed software version numbers.", required=True)
    parser.add_argument("-a", "--alleles", help="<Required> MHC Alleles", required=True, type=str)
    parser.add_argument(
        "-r",
        "--reference",
        help="Reference, retrieved information will be based on this ensembl version",
        required=False,
        default="GRCh37",
        choices=["GRCh37", "GRCh38","COVID"],
    )
    parser.add_argument(
        "-f", "--filter_self", help="Filter peptides against human proteom", required=False, action="store_true"
    )
    parser.add_argument(
        "-wt",
        "--wild_type",
        help="Add wild type sequences of mutated peptides to output",
        required=False,
        action="store_true",
    )
    parser.add_argument(
        "-fo", "--fasta_output", help="Create FASTA file with protein sequences", required=False, action="store_true"
    )
    parser.add_argument("-rp", "--reference_proteome", help="Reference proteome for self-filtering", required=False)
    parser.add_argument("-gr", "--gene_reference", help="List of gene IDs for ID mapping.", required=False)
    parser.add_argument("-pq", "--protein_quantification", help="File with protein quantification values")
    parser.add_argument("-ge", "--gene_expression", help="File with expression analysis results")
    parser.add_argument(
        "-de", "--diff_gene_expression", help="File with differential expression analysis results (DESeq2)"
    )
    parser.add_argument(
        "-li",
        "--ligandomics_id",
        help="Comma separated file with peptide sequence, score and median intensity of a ligandomics identification run.",
    )

    parser.add_argument(
        "-var",
        "--variant_lineage",
        required=False,
        default="Unknown",
        help="Variant lineage",
    )

    parser.add_argument("-v", "--version", help="Script version", action="version", version=VERSION)
    args = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit("Provide at least one argument to epaa.py.")

    filehandler = logging.FileHandler("{}_prediction.log".format(args.identifier))
    filehandler.setLevel(logging.DEBUG)
    filehandler.setFormatter(formatter)
    logger.addHandler(filehandler)

    logger.info("Running Epitope Prediction And Annotation version: " + str(VERSION))
    logger.info("Starting predictions at " + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    metadata = []
    proteins = []
    #references = {"GRCh37": "https://feb2014.archive.ensembl.org", "GRCh38": "http://apr2018.archive.ensembl.org"}
    references = {"GRCh38": "http://jul2023.archive.ensembl.org"}
    global transcriptProteinMap
    global transcriptSwissProtMap

    if args.somatic_mutations.endswith(".vcf"):
        logger.info("Running epaa for variants...")
        variant_list, transcripts, metadata = read_vcf(args.somatic_mutations)

    else:
        logger.error("Please provide a vcf.")

        transcripts = list(set(transcripts))
        transcriptProteinMap, transcriptSwissProtMap = get_protein_ids_for_transcripts(
            ID_SYSTEM_USED, transcripts, references[args.reference], args.reference
        )


    # get the alleles
    # alleles = FileReader.read_lines(args.alleles, in_type=Allele)
    alleles = [Allele(a) for a in args.alleles.split(";")]

    # initialize MartsAdapter, GRCh37 or GRCh38 based
    ma = MartsAdapter(biomart=references[args.reference])

    # create protein db instance for filtering self-peptides
    up_db = UniProtDB("sp")
    if args.filter_self:
        logger.info("Reading human proteome")

        if os.path.isdir(args.reference_proteome):
            for filename in os.listdir(args.reference_proteome):
                if filename.endswith(".fasta") or filename.endswith(".fsa"):
                    up_db.read_seqs(os.path.join(args.reference_proteome, filename))
        else:
            up_db.read_seqs(args.reference_proteome)
        

    selected_methods = [item.split("-")[0] if "mhcnuggets" not in item else item for item in args.tools.split(",")]
    with open(args.versions, "r") as versions_file:
        tool_version = [(row[0].split()[0], str(row[1])) for row in csv.reader(versions_file, delimiter=":")]
        # NOTE this needs to be updated, if a newer version will be available via epytope and should be used in the future
        tool_version.append(("syfpeithi", "1.0"))
        # get for each selected method the corresponding tool version
        methods = {
            method.lower().strip(): version.strip()
            for tool, version in tool_version
            for method in selected_methods
            if tool.lower() in method.lower()
        }

    for method, version in methods.items():
        if version not in EpitopePredictorFactory.available_methods()[method]:
            raise ValueError("The specified version " + version + " for " + method + " is not supported by epytope.")

    thresholds = {
        "syfpeithi": 50,
        "mhcflurry": 500,
        "mhcnuggets-class-1": 500,
        "mhcnuggets-class-2": 500,
        "netmhc": 500,
        "netmhcpan": 500,
        "netmhcii": 500,
        "netmhciipan": 500,
    }
    # Define binders based on the rank metric for netmhc family tools
    # NOTE these recommended thresholds might change in the future with new versions of the tools
    if "netmhc" in "".join(methods.keys()) and not args.use_affinity_thresholds:
        thresholds.update({"netmhc": 2, "netmhcpan": 2, "netmhcii": 10, "netmhciipan": 5})

    if args.tool_thresholds:
        with open(args.tool_thresholds, "r") as json_file:
            threshold_file = json.load(json_file)
            for tool, thresh in threshold_file.items():
                if tool in thresholds.keys():
                    thresholds[tool] = thresh
                else:
                    raise ValueError("Tool " + tool + " in specified threshold file is not supported")
    
    # Distinguish between prediction for peptides and variants
    if args.peptides:
        pred_dataframes, statistics = make_predictions_from_peptides(
            peptides, methods, thresholds, args.use_affinity_thresholds, alleles, up_db, args.identifier, metadata
        )
    else:
        pred_dataframes, statistics, all_peptides_filtered, proteins = make_predictions_from_variants(
            variant_list,
            methods,
            thresholds,
            args.use_affinity_thresholds,
            alleles,
            int(args.min_length),
            int(args.max_length) + 1,
            ma,
            up_db,
            args.identifier,
            metadata,
            transcriptProteinMap,
        )
    # concat dataframes for all peptide lengths
    try:
        complete_df = pd.concat(pred_dataframes, sort=True)
        # replace method names with method names with version
        complete_df["method"] = complete_df["method"].apply(lambda x: x.lower() + "-" + methods[x.lower()])
        complete_df["variant_lineage"] = args.variant_lineage
        predictions_available = True
    except:
        complete_df = pd.DataFrame()
        predictions_available = False
        logger.error("No predictions available.")

    # include wild type sequences to dataframe if specified
    if args.wild_type:
        if args.peptides:
            logger.warning("Wildtype sequence generation not available with peptide input.")
            pass
        wt_sequences = generate_wt_seqs(all_peptides_filtered)
        complete_df["wt sequence"] = complete_df.apply(
            lambda row: create_wt_seq_column_value(row, wt_sequences), axis=1
        )
        columns_tiles = [
            "sequence",
            "refseq",
            "wt sequence",
            "length",
            "chr",
            "pos",
            "gene",
            "transcripts",
            "proteins",
            "variant type",
            "method",
        ]
    # Change the order (the index) of the columns
    else:
        columns_tiles = [
            "sequence",
            "refseq",
            "length",
            "chr",
            "pos",
            "gene",
            "transcripts",
            "proteins",
            "variant type",
            "method",
        ]

    for c in complete_df.columns:
        if c not in columns_tiles:
            columns_tiles.append(c)
            
    complete_df = complete_df.reindex(columns=columns_tiles)

    #TODO checar se nao ta pegando de refseq
    #TODO fazer para binder_changes
    binder_cols = [col for col in complete_df.columns if ("binder" in col) and ("change" not in col) and ("refseq" not in col)]
    binder_changes_cols = [col for col in complete_df.columns if ("binder_changes" in col)]


    binders = []
    non_binders = []
    pos_predictions = []
    neg_predictions = []

    change_binders = []
    

    for i, r in complete_df.iterrows():
        binder = False
        change = False
        for c in binder_cols:
            if r[c] is True:
                binder = True
                continue
        if binder:
            binders.append(str(r["sequence"]))
            pos_predictions.append(str(r["sequence"]))
        else:
            neg_predictions.append(str(r["sequence"]))
            if str(r["sequence"]) not in binders:
                non_binders.append(str(r["sequence"]))

        for c in binder_changes_cols:
            if r[c] is True:
                change = True
                continue
        if change:
            change_binders .append(str(r["sequence"]))
        
    # parse protein quantification results, annotate proteins for samples
    if args.protein_quantification is not None:
        protein_quant = read_protein_quant(args.protein_quantification)
        first_entry = protein_quant[protein_quant.keys()[0]]
        for k in first_entry.keys():
            complete_df["{} log2 protein LFQ intensity".format(k)] = complete_df.apply(
                lambda row: create_quant_column_value_for_result(row, protein_quant, transcriptSwissProtMap, k), axis=1
            )
    # parse (differential) expression analysis results, annotate features (genes/transcripts)
    if args.gene_expression is not None:
        fold_changes = read_diff_expression_values(args.gene_expression)
        gene_id_lengths = {}
        col_name = "RNA expression (RPKM)"

        with open(args.gene_reference, "r") as gene_list:
            for l in gene_list:
                ids = l.split("\t")
                gene_id_in_df = complete_df.iloc[1]["gene"]
                if "ENSG" in gene_id_in_df:
                    gene_id_lengths[ids[0]] = float(ids[2].strip())
                else:
                    gene_id_lengths[ids[1]] = float(ids[2].strip())
        deseq = False
        # add column to result dataframe
        complete_df[col_name] = complete_df.apply(
            lambda row: create_expression_column_value_for_result(row, fold_changes, deseq, gene_id_lengths), axis=1
        )
    if args.diff_gene_expression is not None:
        gene_id_lengths = {}
        fold_changes = read_diff_expression_values(args.diff_gene_expression)
        col_name = "RNA normal_vs_tumor.log2FoldChange"
        deseq = True

        # add column to result dataframe
        complete_df[col_name] = complete_df.apply(
            lambda row: create_expression_column_value_for_result(row, fold_changes, deseq, gene_id_lengths), axis=1
        )
    # parse ligandomics identification results, annotate peptides for samples
    if args.ligandomics_id is not None:
        lig_id = read_lig_ID_values(args.ligandomics_id)
        # add columns to result dataframe
        complete_df["ligand score"] = complete_df.apply(
            lambda row: create_ligandomics_column_value_for_result(row, lig_id, 0, False), axis=1
        )
        complete_df["ligand intensity"] = complete_df.apply(
            lambda row: create_ligandomics_column_value_for_result(row, lig_id, 1, False), axis=1
        )

        if args.wild_type != None:
            complete_df["wt ligand score"] = complete_df.apply(
                lambda row: create_ligandomics_column_value_for_result(row, lig_id, 0, True), axis=1
            )
            complete_df["wt ligand intensity"] = complete_df.apply(
                lambda row: create_ligandomics_column_value_for_result(row, lig_id, 1, True), axis=1
            )
    # write mutated protein sequences to fasta file
    if args.fasta_output and predictions_available:
        with open("{}_prediction_proteins.fasta".format(args.identifier), "w") as protein_outfile:
            for p in proteins:
                variants = []
                for v in p.vars:
                    variants = variants + p.vars[v]
                c = [x.coding.values() for x in variants]
                cf = list(itertools.chain.from_iterable(c))
                cds = ",".join([y.cdsMutationSyntax for y in set(cf)])
                aas = ",".join([y.aaMutationSyntax for y in set(cf)])
                protein_outfile.write(">{}:{}:{}\n".format(p.transcript_id, aas, cds))
                protein_outfile.write("{}\n".format(str(p)))

    complete_df["binder"] = complete_df[[col for col in complete_df.columns if "binder" in col]].any(axis=1)

    # write dataframe to tsv
    complete_df.fillna("")
    if predictions_available:
        complete_df.to_csv("{}_prediction_result.tsv".format(args.identifier), "\t", index=False)

    statistics["tool_thresholds"] = thresholds
    statistics["number_of_predictions"] = len(complete_df)
    statistics["number_of_binders"] = len(pos_predictions)
    statistics["number_of_nonbinders"] = len(neg_predictions)
    statistics["number_of_changed"] = len(change_binders)
    statistics["number_of_unique_binders"] = list(set(binders))
    statistics["number_of_unique_nonbinders"] = list(set(non_binders) - set(binders))
    statistics["number_of_unique_changed_binders"] = list(set(change_binders))

    with open("{}_report.json".format(args.identifier), "w") as json_out:
        json.dump(statistics, json_out)

    logger.info("Finished predictions at " + str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))


if __name__ == "__main__":
    __main__()
