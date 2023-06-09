# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import shutil

import pandas as pd
import qiime2
import skbio
from q2_types.feature_data import DNAFASTAFormat, DNAIterator, ProteinFASTAFormat
from q2_types.feature_data._transformer import ProteinIterator
from q2_types_genomics.genome_data import GenesDirectoryFormat, ProteinsDirectoryFormat
from skbio import DNA, Protein

from q2_amr.types import CARDAnnotationDirectoryFormat

from ..plugin_setup import plugin
from ._format import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationJSONFormat,
    CARDAnnotationTXTFormat,
    CARDDatabaseFormat,
    CARDGeneAnnotationDirectoryFormat,
)


@plugin.register_transformer
def _1(data: CARDDatabaseFormat) -> pd.DataFrame:
    ff = pd.read_json(str(data)).transpose()
    return ff


@plugin.register_transformer
def _2(data: pd.DataFrame) -> CARDDatabaseFormat:
    ff = CARDDatabaseFormat()
    with ff.open() as fh:
        data.transpose().to_json(fh)
    return ff


@plugin.register_transformer
def _3(data: CARDDatabaseFormat) -> DNAFASTAFormat:
    path = str(data)
    file = _read_from_card_file(path, "dna")
    return file


@plugin.register_transformer
def _4(data: CARDDatabaseFormat) -> ProteinFASTAFormat:
    path = str(data)
    file = _read_from_card_file(path, "protein")
    return file


@plugin.register_transformer
def _5(data: CARDDatabaseFormat) -> DNAIterator:
    path = str(data)
    generator = _read_from_card_generator(path, "dna")
    return DNAIterator(generator)


@plugin.register_transformer
def _6(data: CARDDatabaseFormat) -> ProteinIterator:
    path = str(data)
    generator = _read_from_card_generator(path, "protein")
    return ProteinIterator(generator)


def _read_from_card_file(path, seq_type):
    with open(str(path), "rb") as f:
        db = json.load(f)
    genomes = ProteinFASTAFormat() if seq_type == "protein" else DNAFASTAFormat()
    with open(str(genomes), "a") as fin:
        for key1, value in db.items():
            if "model_sequences" not in value:
                continue
            key2 = list(value["model_sequences"]["sequence"].keys())[0]
            sequence = extract_sequence(seq_type, key1, key2, db)
            if sequence:
                skbio.io.write(sequence, format="fasta", into=fin)
            else:
                continue
        return genomes


def _read_from_card_generator(path, seq_type):
    with open(str(path), "rb") as f:
        db = json.load(f)
    for key1, value in db.items():
        if "model_sequences" not in value:
            continue
        key2 = list(value["model_sequences"]["sequence"].keys())[0]
        sequence = extract_sequence(seq_type, key1, key2, db)
        if sequence:
            yield sequence
        else:
            continue


def extract_sequence(seq_type, key1, key2, db):
    if seq_type == "protein":
        obj = Protein
    else:
        obj = DNA

    sequence = db[key1]["model_sequences"]["sequence"][key2][f"{seq_type}_sequence"][
        "sequence"
    ]
    if sequence is None:
        return None
    sequence_object = obj(
        db[key1]["model_sequences"]["sequence"][key2][f"{seq_type}_sequence"][
            "sequence"
        ]
    )
    aro = db[key1]["ARO_accession"]
    accession = db[key1]["model_sequences"]["sequence"][key2][f"{seq_type}_sequence"][
        "accession"
    ]
    model_name = db[key1]["model_name"]
    ncbi_taxonomy = db[key1]["model_sequences"]["sequence"][key2]["NCBI_taxonomy"][
        "NCBI_taxonomy_name"
    ]
    sequence_object.metadata["description"] = f"[{ncbi_taxonomy}]"
    if seq_type == "protein":
        sequence_object.metadata["id"] = f"gb|{accession}|ARO:{aro}|{model_name}"
    else:
        strand = db[key1]["model_sequences"]["sequence"][key2]["dna_sequence"]["strand"]
        fmin = db[key1]["model_sequences"]["sequence"][key2]["dna_sequence"]["fmin"]
        fmax = db[key1]["model_sequences"]["sequence"][key2]["dna_sequence"]["fmax"]
        sequence_object.metadata[
            "id"
        ] = f"gb|{accession}|{strand}|{fmin}-{fmax}|ARO:{aro}|{model_name}"
    return sequence_object


@plugin.register_transformer
def _7(data: pd.DataFrame) -> CARDAnnotationTXTFormat:
    ff = CARDAnnotationTXTFormat()
    with ff.open() as fh:
        data.to_csv(fh, index=False, sep="\t")
    return ff


@plugin.register_transformer
def _8(data: CARDAnnotationTXTFormat) -> pd.DataFrame:
    file = pd.read_csv(str(data), sep="\t")
    return file


@plugin.register_transformer
def _9(data: dict) -> CARDAnnotationJSONFormat:
    ff = CARDAnnotationJSONFormat()
    with ff.open() as fh:
        json.dump(data, fh)
    return ff


@plugin.register_transformer
def _10(data: CARDAnnotationDirectoryFormat) -> GenesDirectoryFormat:
    genes_directory = GenesDirectoryFormat()
    create_dir_structure(data, "DNA", genes_directory)
    return genes_directory


@plugin.register_transformer
def _11(data: CARDAnnotationDirectoryFormat) -> ProteinsDirectoryFormat:
    proteins_directory = ProteinsDirectoryFormat()
    create_dir_structure(data, "Protein", proteins_directory)
    return proteins_directory


def create_dir_structure(data, seq_type, genes_protein_directory):
    annotation_dir = str(data)
    for sample in os.listdir(annotation_dir):
        for bin in os.listdir(os.path.join(annotation_dir, sample)):
            for file in os.listdir(os.path.join(annotation_dir, sample, bin)):
                if file.endswith(".txt"):
                    txt_file_path = os.path.join(annotation_dir, sample, bin, file)
                    os.makedirs(
                        os.path.join(str(genes_protein_directory), sample),
                        exist_ok=True,
                    )
                    fasta = card_annotation_df_to_fasta(txt_file_path, seq_type)
                    filename = (
                        f"{bin}_genes.fasta"
                        if seq_type == "DNA"
                        else f"{bin}_proteins.fasta"
                    )
                    shutil.move(
                        str(fasta),
                        os.path.join(str(genes_protein_directory), sample, filename),
                    )


def card_annotation_df_to_fasta(txt_file_path: str, seq_type: str):
    annotation_df = pd.read_csv(txt_file_path, sep="\t")
    fasta_format, sequence_class = (
        (DNAFASTAFormat(), DNA)
        if seq_type == "DNA"
        else (ProteinFASTAFormat(), Protein)
    )
    with open(str(fasta_format), "a") as fasta_file:
        for index, row in annotation_df.iterrows():
            sequence_object = sequence_class(row[f"Predicted_{seq_type}"])
            sequence_object.metadata["id"] = row["ORF_ID"]
            sequence_object.metadata["description"] = row["ARO"]
            skbio.io.write(sequence_object, format="fasta", into=fasta_file)
    return fasta_format


def read_mapping_data(data_path, variant):
    df_list = []
    for samp in os.listdir(str(data_path)):
        file_path = os.path.join(
            str(data_path), samp, f"{samp}.{variant}_mapping_data.txt"
        )
        df = pd.read_csv(file_path, sep="\t")
        df.insert(0, "Sample Name", samp)
        df_list.append(df)
    mapping_data_cat = pd.concat(df_list, axis=0)
    mapping_data_cat.reset_index(inplace=True, drop=True)
    mapping_data_cat.index.name = "id"
    mapping_data_cat.index = mapping_data_cat.index.astype(str)
    metadata = qiime2.Metadata(mapping_data_cat)
    return metadata


@plugin.register_transformer
def _12(data: CARDAlleleAnnotationDirectoryFormat) -> qiime2.Metadata:
    return read_mapping_data(data, "allele")


@plugin.register_transformer
def _13(data: CARDGeneAnnotationDirectoryFormat) -> qiime2.Metadata:
    return read_mapping_data(data, "gene")
