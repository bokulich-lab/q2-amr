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
import skbio
from q2_types_genomics.genome_data import GenesDirectoryFormat

from q2_amr.types import CARDAnnotationDirectoryFormat
from q2_types.feature_data import (
    DNAFASTAFormat, DNAIterator, ProteinFASTAFormat)
from q2_types.feature_data._transformer import ProteinIterator
from skbio import Protein, DNA

from ..plugin_setup import plugin
from ._format import CARDDatabaseFormat, CARDAnnotationTXTFormat, CARDAnnotationJSONFormat


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
    file = _read_from_card_file(path, 'dna')
    return file

@plugin.register_transformer
def _12(data: CARDDatabaseFormat) -> ProteinFASTAFormat:
    path = str(data)
    file = _read_from_card_file(path, 'protein')
    return file

@plugin.register_transformer
def _13(data: CARDDatabaseFormat) -> DNAIterator:
    path = str(data)
    generator = _read_from_card_generator(path, 'dna')
    return DNAIterator(generator)

@plugin.register_transformer
def _14(data: CARDDatabaseFormat) -> ProteinIterator:
    path = str(data)
    generator = _read_from_card_generator(path, 'protein')
    return ProteinIterator(generator)


def _read_from_card_file(path, seq_type):
    with open(str(path), 'rb') as f:
        db = json.load(f)
    genomes = ProteinFASTAFormat() if seq_type == 'protein' else DNAFASTAFormat()
    with open(str(genomes), 'a') as fin:
        for key1, value in db.items():
            if 'model_sequences' not in value:
                continue
            key2 = list(value['model_sequences']['sequence'].keys())[0]
            sequence = extract_sequence(seq_type, key1, key2, db)
            if sequence:
                skbio.io.write(sequence, format='fasta', into=fin)
            else:
                continue
        return genomes


def _read_from_card_generator(path, seq_type):
    with open(str(path), 'rb') as f:
        db = json.load(f)
    for key1, value in db.items():
        if 'model_sequences' not in value:
            continue
        key2 = list(value['model_sequences']['sequence'].keys())[0]
        sequence = extract_sequence(seq_type, key1, key2, db)
        if sequence:
            yield sequence
        else:
            continue

def extract_sequence(seq_type, key1, key2, db):
    if seq_type == 'protein':
        obj = Protein
    else:
        obj = DNA

    sequence = db[key1]['model_sequences']['sequence'][key2][f'{seq_type}_sequence']['sequence']
    if sequence is None:
        return None
    sequence_object = obj(db[key1]['model_sequences']['sequence'][key2][f'{seq_type}_sequence']['sequence'])
    aro = db[key1]['ARO_accession']
    accession = db[key1]['model_sequences']['sequence'][key2][f'{seq_type}_sequence']['accession']
    model_name = db[key1]['model_name']
    ncbi_taxonomy = db[key1]['model_sequences']['sequence'][key2]['NCBI_taxonomy']['NCBI_taxonomy_name']
    sequence_object.metadata['description'] = f'[{ncbi_taxonomy}]'
    if seq_type == 'protein':
        sequence_object.metadata['id'] = f'gb|{accession}|ARO:{aro}|{model_name}'
    else:
        strand = db[key1]['model_sequences']['sequence'][key2]['dna_sequence']['strand']
        fmin = db[key1]['model_sequences']['sequence'][key2]['dna_sequence']['fmin']
        fmax = db[key1]['model_sequences']['sequence'][key2]['dna_sequence']['fmax']
        sequence_object.metadata['id'] = f'gb|{accession}|{strand}|{fmin}-{fmax}|ARO:{aro}|{model_name}'
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
def _16(data: dict) -> CARDAnnotationJSONFormat:
    ff = CARDAnnotationJSONFormat()
    with ff.open() as fh:
        json.dump(data, fh)
    return ff

@plugin.register_transformer
def _16(data: CARDAnnotationDirectoryFormat) -> GenesDirectoryFormat:
    genes_directory = GenesDirectoryFormat()
    annotation_dir = str(data)
    for sample in os.listdir(annotation_dir):
        for bin in os.listdir(os.path.join(annotation_dir, sample)):
            for file in os.listdir(os.path.join(annotation_dir, sample, bin)):
                if file.endswith('.txt'):
                    txt_file_path = os.path.join(annotation_dir, sample, bin, file)
                    os.makedirs(os.path.join(str(genes_directory), sample), exist_ok=True)
                    fasta = card_annotation_df_to_fasta(txt_file_path, 'gene')
                    shutil.move(str(fasta), os.path.join(str(genes_directory), sample, f'{bin}_genes.fasta'))
    return genes_directory


def card_annotation_df_to_fasta(txt_file_path: str, sort):
    annotation_df = pd.read_csv(txt_file_path, sep='\t')
    if sort is "gene":
        fasta = DNAFASTAFormat()
        with open(str(fasta), 'a') as dnaf:
            for index, row in annotation_df.iterrows():
                dna_object = DNA(row['Predicted_DNA'])
                dna_object.metadata['id'] = row['ORF_ID']
                dna_object.metadata['description'] = row['ARO']
                skbio.io.write(dna_object, format='fasta', into=dnaf)
    else:
        fasta = ProteinFASTAFormat()
        with open(str(fasta), 'a') as proteinf:
            for index, row in annotation_df.iterrows():
                protein_object = Protein(row['Predicted_Protein'])
                protein_object.metadata['id'] = row['ORF_ID']
                protein_object.metadata['description'] = row['ARO']
                skbio.io.write(protein_object, format='fasta', into=proteinf)
    return fasta
