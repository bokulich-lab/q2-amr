# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os

import pandas as pd
import qiime2
import skbio
from q2_types.feature_data import DNAFASTAFormat, DNAIterator, ProteinFASTAFormat
from q2_types.feature_data._transformer import ProteinIterator
from skbio import DNA, Protein

from ..plugin_setup import plugin
from ._format import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationJSONFormat,
    CARDAnnotationTXTFormat,
    CARDDatabaseFormat,
    CARDGeneAnnotationDirectoryFormat,
)


@plugin.register_transformer
def _9(data: CARDDatabaseFormat) -> pd.DataFrame:
    ff = pd.read_json(str(data)).transpose()
    return ff


@plugin.register_transformer
def _10(data: pd.DataFrame) -> CARDDatabaseFormat:
    ff = CARDDatabaseFormat()
    with ff.open() as fh:
        data.transpose().to_json(fh)
    return ff


@plugin.register_transformer
def _11(data: CARDDatabaseFormat) -> DNAFASTAFormat:
    path = str(data)
    file = _read_from_card_file(path, "dna")
    return file


@plugin.register_transformer
def _12(data: CARDDatabaseFormat) -> ProteinFASTAFormat:
    path = str(data)
    file = _read_from_card_file(path, "protein")
    return file


@plugin.register_transformer
def _13(data: CARDDatabaseFormat) -> DNAIterator:
    path = str(data)
    generator = _read_from_card_generator(path, "dna")
    return DNAIterator(generator)


@plugin.register_transformer
def _14(data: CARDDatabaseFormat) -> ProteinIterator:
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
def _15(data: pd.DataFrame) -> CARDAnnotationTXTFormat:
    ff = CARDAnnotationTXTFormat()
    with ff.open() as fh:
        data.to_csv(fh, index=False, sep="\t")
    return ff


@plugin.register_transformer
def _16(data: CARDAnnotationTXTFormat) -> pd.DataFrame:
    file = pd.read_csv(str(data), sep="\t")
    return file


@plugin.register_transformer
def _17(data: dict) -> CARDAnnotationJSONFormat:
    ff = CARDAnnotationJSONFormat()
    with ff.open() as fh:
        json.dump(data, fh)
    return ff


@plugin.register_transformer
def _18(data: CARDAlleleAnnotationDirectoryFormat) -> qiime2.Metadata:
    df_list = []
    for samp in os.listdir(str(data)):
        df = pd.read_csv(
            os.path.join(str(data), samp, f"{samp}.allele_mapping_data.txt"), sep="\t"
        )
        df.insert(0, "Sample Name", samp)
        df_list.append(df)
    mapping_data_cat = pd.concat(df_list, axis=0)
    mapping_data_cat.reset_index(inplace=True, drop=True)
    mapping_data_cat.index.name = "id"
    mapping_data_cat.index = mapping_data_cat.index.astype(str)
    metadata = qiime2.Metadata(mapping_data_cat)
    return metadata


@plugin.register_transformer
def _19(data: CARDGeneAnnotationDirectoryFormat) -> qiime2.Metadata:
    df_list = []
    for samp in os.listdir(str(data)):
        df = pd.read_csv(
            os.path.join(str(data), samp, f"{samp}.gene_mapping_data.txt"), sep="\t"
        )
        df.insert(0, "Sample Name", samp)
        df_list.append(df)
    mapping_data_cat = pd.concat(df_list, axis=0)
    mapping_data_cat.reset_index(inplace=True, drop=True)
    mapping_data_cat.index.name = "id"
    mapping_data_cat.index = mapping_data_cat.index.astype(str)
    metadata = qiime2.Metadata(mapping_data_cat)
    return metadata
