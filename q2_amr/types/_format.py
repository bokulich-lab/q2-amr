# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
import qiime2.plugin.model as model
from qiime2.plugin import ValidationError


class CARDDatabaseFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        HEADER = ['model_id', 'model_name', 'model_type', 'model_type_id', 'model_description', 'model_param', 'model_sequences', 'ARO_accession', 'ARO_id', 'ARO_name', 'CARD_short_name', 'ARO_description', 'ARO_category', 'description', 'access']
        HEADER2 = ['model_id', 'model_name', 'model_type', 'model_type_id', 'model_description', 'model_param', 'model_sequences', 'ARO_accession', 'ARO_id', 'ARO_name', 'ARO_description', 'ARO_category', 'description', 'access']
        card_df = pd.read_json(str(self)).transpose()
        header = list(card_df.columns)
        if header != HEADER and header != HEADER2:
            raise ValidationError(
                "Header line does not match CARDDatabase format. Must consist of "
                "the following values: " + ', '.join(HEADER) +
                ".\n\nFound instead: " + ', '.join(header))

    def _validate_(self, level):
        self._validate()


CARDDatabaseDirectoryFormat = model.SingleFileDirectoryFormat(
    'CARDDatabaseDirectoryFormat', 'card.json',
    CARDDatabaseFormat)


class CARDAnnotationtxtFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        HEADER = ['ORF_ID', 'Contig', 'Start', 'Stop', 'Orientation', 'Cut_Off', 'Pass_Bitscore', 'Best_Hit_Bitscore', 'Best_Hit_ARO', 'Best_Identities', 'ARO', 'Model_type', 'SNPs_in_Best_Hit_ARO', 'Other_SNPs', 'Drug Class', 'Resistance Mechanism', 'AMR Gene Family', 'Predicted_DNA', 'Predicted_Protein', 'CARD_Protein_Sequence', 'Percentage Length of Reference Sequence', 'ID', 'Model_ID', 'Nudged', 'Note']
        df = pd.read_csv(str(self), sep="\t")

        header = list(df.columns)
        if header != HEADER:
            raise ValidationError(
                "Header line does not match CARDAnnotation format. Must consist of "
                "the following values: " + ', '.join(HEADER) +
                ".\n\nFound instead: " + ', '.join(header))

    def _validate_(self, level):
        self._validate()


CARDAnnotationtxtDirectoryFormat = model.SingleFileDirectoryFormat(
    'CARDAnnotationDirectoryFormat', 'amr_annotation.txt', CARDAnnotationtxtFormat)


class CARDAnnotationjsonFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        x = 1
        if x is 2:
            raise ValidationError(
                "Header line does not match CARDAnnotation format. Must consist of "
                "the following values: " + ', '.join('n') +
                ".\n\nFound instead: " + ', '.join('n'))

    def _validate_(self, level):
        self._validate()


CARDAnnotationjsonDirectoryFormat = model.SingleFileDirectoryFormat(
    'CARDAnnotationjsonDirectoryFormat', 'amr_annotation.json', CARDAnnotationjsonFormat)
