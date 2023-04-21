# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
from copy import copy

import pandas as pd
import qiime2.plugin.model as model
from qiime2.plugin import ValidationError


class CARDDatabaseFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        header_exp = ['model_id', 'model_name', 'model_type', 'model_type_id', 'model_description', 'model_param', 'model_sequences', 'ARO_accession', 'ARO_id', 'ARO_name', 'CARD_short_name', 'ARO_description', 'ARO_category', 'description', 'access']

        header_exp_2 = copy(header_exp)
        header_exp_2.pop(10)
        card_df = pd.read_json(str(self)).transpose()
        header_obs = list(card_df.columns)
        if header_obs != header_exp and header_obs != header_exp_2:
            raise ValidationError(
                "Header line does not match CARDDatabase format. Must consist of "
                "the following values: " + ', '.join(header_exp) +
                ".\n\nFound instead: " + ', '.join(header_obs))

    def _validate_(self, level):
        self._validate()


CARDDatabaseDirectoryFormat = model.SingleFileDirectoryFormat(
    'CARDDatabaseDirectoryFormat', 'card.json',
    CARDDatabaseFormat)


class CARDAnnotationTXTFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        header_exp = ['ORF_ID', 'Contig', 'Start', 'Stop', 'Orientation', 'Cut_Off', 'Pass_Bitscore', 'Best_Hit_Bitscore',
                  'Best_Hit_ARO', 'Best_Identities', 'ARO', 'Model_type', 'SNPs_in_Best_Hit_ARO', 'Other_SNPs',
                  'Drug Class', 'Resistance Mechanism', 'AMR Gene Family', 'Predicted_DNA', 'Predicted_Protein',
                  'CARD_Protein_Sequence', 'Percentage Length of Reference Sequence', 'ID', 'Model_ID', 'Nudged',
                  'Note']
        df = pd.read_csv(str(self), sep="\t")

        header_obs = list(df.columns)
        if header_obs != header_exp:
            raise ValidationError(
                "Header line does not match CARDAnnotation format. Must consist of "
                "the following values: " + ', '.join(header_exp) +
                ".\n\nFound instead: " + ', '.join(header_obs))

    def _validate_(self, level):
        self._validate()


class CARDAnnotationJSONFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        keys_exp = ['bit_score', 'sequence_from_broadstreet', 'query_start', 'ARO_name', 'sequence_from_db', 'max_identities', 'orf_prot_sequence', 'query_snp', 'match', 'orf_strand', 'model_id', 'evalue', 'model_name', 'model_type', 'model_type_id', 'orf_from', 'ARO_accession', 'partial', 'orf_end', 'dna_sequence_from_broadstreet', 'perc_identity', 'query_end', 'cvterm_id', 'type_match', 'pass_evalue', 'pass_bitscore', 'ARO_category', 'snp', 'orf_dna_sequence', 'query', 'orf_start']

        with open(str(self), 'r') as f:
            json_str = f.read()
            json_data = json.loads(json_str)

        # Traverse the nested data structure to access the keys
        keys_obs = []
        for k, v in json_data.items():
            for sub_k, sub_v in v.items():
                keys_obs.extend(sub_v.keys())
        # remove duplicates
        keys_obs = set(keys_obs)

        if keys_obs != set(keys_exp):
            raise ValidationError(
                "Dict keys do not match CARDAnnotation format. Must consist of "
                "the following values: " + ', '.join(keys_exp) +
                ".\n\nFound instead: " + ', '.join(keys_obs))

    def _validate_(self, level):
        self._validate()

class CARDAnnotationDirectoryFormat(model.DirectoryFormat):
    json = model.File(r'amr_annotation.json', format=CARDAnnotationJSONFormat)
    txt = model.File(r'amr_annotation.txt', format=CARDAnnotationTXTFormat)
