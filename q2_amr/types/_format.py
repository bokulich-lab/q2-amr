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
from q2_types_genomics.per_sample_data._format import MultiDirValidationMixin
from qiime2.plugin import ValidationError


class CARDDatabaseFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        header_exp = ['model_id', 'model_name', 'model_type', 'model_type_id', 'model_description', 'model_param',
                      'model_sequences', 'ARO_accession', 'ARO_id', 'ARO_name', 'CARD_short_name', 'ARO_description',
                      'ARO_category', 'description', 'access']

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
        header_exp = ['ORF_ID', 'Contig', 'Start', 'Stop', 'Orientation', 'Cut_Off', 'Pass_Bitscore',
                      'Best_Hit_Bitscore',
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
        keys_exp = ['match', 'cvterm_id', 'orf_prot_sequence', 'model_id', 'ARO_category', 'orf_start', 'ARO_accession',
                    'evalue', 'sequence_from_broadstreet', 'query', 'model_type_id', 'model_type', 'bit_score',
                    'sequence_from_db', 'query_end', 'orf_dna_sequence', 'pass_bitscore', 'orf_end', 'pass_evalue',
                    'query_start', 'perc_identity', 'type_match', 'max_identities', 'orf_from', 'ARO_name',
                    'model_name', 'orf_strand']
        keys_obs = []
        with open(str(self), 'r') as f:
            json_str = f.read()
            json_data = json.loads(json_str)

        for k, v in json_data.items():
            for sub_k, sub_v in v.items():
                keys_obs.extend(sub_v.keys())

        if keys_obs and not set(keys_exp).issubset(set(keys_obs)):
            raise ValidationError(
                "Dict keys do not match CARDAnnotation format. Must consist of "
                "the following values: " + ', '.join(keys_exp) +
                ".\n\nFound instead: " + ', '.join(set(keys_obs)))

    def _validate_(self, level):
        self._validate()


class CARDAnnotationDirectoryFormat(MultiDirValidationMixin, model.DirectoryFormat):
    json = model.FileCollection(r'.+amr_annotation.json$', format=CARDAnnotationJSONFormat)
    txt = model.FileCollection(r'.+amr_annotation.txt$', format=CARDAnnotationTXTFormat)

    @json.set_path_maker
    def json_path_maker(self, sample_id, bin_id):
        return f"{sample_id}/{bin_id}/{sample_id}_{bin_id}_amr_annotation.json"

    @txt.set_path_maker
    def txt_path_maker(self, sample_id, bin_id):
        return f"{sample_id}/{bin_id}/{sample_id}_{bin_id}_amr_annotation.txt"
