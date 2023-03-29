# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import tarfile

import pandas as pd
import pandas.util.testing as pdt
import numpy as np
import tempfile
import pkg_resources
import requests
from pandas._testing import assert_frame_equal

from qiime2.plugin import ValidationError
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import (
    FeatureData, DNAIterator, AlignedDNAFASTAFormat, AlignedDNAIterator,
    DNAFASTAFormat, ProteinIterator, ProteinFASTAFormat)  #RNAFASTAFormat
from skbio import DNA, Protein
from skbio.sequence import GrammaredSequence

from rescript.card import fetch_card_data, card_annotation, card_annotation_txt_to_fasta

from rescript.types._format import CARDDatabaseFormat, CARDAnnotationFormat

from rescript._utilities import _read_fasta, _read_dna_alignment_fasta
from rescript.types import (SILVATaxonomyFormat, SILVATaxonomyDirectoryFormat,
                            SILVATaxidMapFormat, SILVATaxidMapDirectoryFormat,
                            SILVATaxonomy, SILVATaxidMap)
from rescript.types._transformer import extract_sequence, _read_from_card_file

from unittest.mock import patch, MagicMock, Mock
from Bio import SeqIO

class AMRTypesTestPluginBase(TestPluginBase):
    package = 'rescript.types.tests'

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.TemporaryDirectory(
            prefix='rescript-test-temp-')

    def tearDown(self):
        self.temp_dir.cleanup()

    def get_data_path(self, filename):
        return pkg_resources.resource_filename(self.package,
                                               'data/%s' % filename)


class TestCARDDatabaseTypesAndFormats(AMRTypesTestPluginBase):

    def test_card_database_format_validate_positive(self):
        filepath = self.get_data_path('card.json')
        format = CARDDatabaseFormat(filepath, mode='r')
        # These should both just succeed
        format.validate()

    def test_dataframe_to_card_format_transformer(self):
        filepath = self.get_data_path('card.json')
        transformer = self.get_transformer(pd.DataFrame, CARDDatabaseFormat)
        card_df = pd.read_json(filepath).transpose()
        obs = transformer(card_df)
        self.assertIsInstance(obs, CARDDatabaseFormat)

    def test_card_format_to_dataframe_transformer(self):
        filepath = self.get_data_path('card.json')
        transformer = self.get_transformer(CARDDatabaseFormat, pd.DataFrame)
        card_db = CARDDatabaseFormat(filepath, mode='r')
        obs = transformer(card_db)
        self.assertIsInstance(obs, pd.DataFrame)

    def test_protein_card_iterator_transformer(self):
        filepath = self.get_data_path('card.json')
        transformer = self.get_transformer(CARDDatabaseFormat, ProteinIterator)
        card_db = CARDDatabaseFormat(filepath, mode='r')
        obs = transformer(card_db)
        self.assertIsInstance(obs, ProteinIterator)


    def test_dna_card_iterator_transformer(self):
        filepath = self.get_data_path('card.json')
        transformer = self.get_transformer(CARDDatabaseFormat, DNAIterator)
        card_db = CARDDatabaseFormat(filepath, mode='r')
        obs = transformer(card_db)
        test_fasta = DNAFASTAFormat(self.get_data_path('card_test_data/card_test_dna.fasta'), mode='r')
        exp = test_fasta.view(DNAIterator)
        self.assertIsInstance(obs, DNAIterator)
        for e, o in zip(exp, obs):
            self.assertEqual(e, o)



    def test_protein_card_fasta_transformer(self):
        filepath = self.get_data_path('card_test_data/card_test.json')
        transformer = self.get_transformer(CARDDatabaseFormat, ProteinFASTAFormat)
        card_db = CARDDatabaseFormat(filepath, mode='r')
        obs = transformer(card_db)
        obs_dict = {rec.id: rec.seq for rec in SeqIO.parse(str(obs), "fasta")}
        exp_dict = {rec.id: rec.seq for rec in SeqIO.parse('/Users/vinzent/Desktop/bokulich_project/Code/RESCRIPt/rescript/types/tests/data/card_test_data/card_test_protein.fasta',"fasta")}
        self.assertDictEqual(obs_dict, exp_dict)
        self.assertIsInstance(obs, ProteinFASTAFormat)

    def test_dna_card_fasta_transformer(self):
        filepath = self.get_data_path('card_test_data/card_test.json')
        transformer = self.get_transformer(CARDDatabaseFormat, DNAFASTAFormat)
        card_db = CARDDatabaseFormat(filepath, mode='r')
        obs = transformer(card_db)
        obs_dict = {rec.id: rec.seq for rec in SeqIO.parse(str(obs), "fasta")}
        exp_dict = {rec.id: rec.seq for rec in SeqIO.parse('/Users/vinzent/Desktop/bokulich_project/Code/RESCRIPt/rescript/types/tests/data/card_test_data/card_test_dna.fasta', "fasta")}
        self.assertDictEqual(obs_dict, exp_dict)
        self.assertIsInstance(obs, DNAFASTAFormat)


    @patch('requests.get')
    def test_fetch_card_data(self, mock_requests):
        f = open(self.get_data_path('card.tar.bz2'), 'rb')
        mock_response = MagicMock(raw=f)
        mock_requests.return_value = mock_response
        obs = fetch_card_data(version='3.2.6')
        self.assertIsInstance(obs, pd.DataFrame)
        mock_requests.assert_called_once_with('https://card.mcmaster.ca/download/0/broadstreet-v3.2.6.tar.bz2', stream=True)
        exp = pd.read_json(self.get_data_path('card_test_data/card_test.json')).transpose()
        assert_frame_equal(obs, exp)

    @patch('requests.get', side_effect=requests.ConnectionError)
    def test_fetch_card_data_connection_error(self, mock_requests):
        with self.assertRaisesRegex(requests.ConnectionError, 'Network connectivity problems.'):
            fetch_card_data(version='3.2.6')

    @patch('tarfile.open', side_effect=tarfile.ReadError)
    def test_fetch_card_data_tarfile_read_error(self, mock_requests):
        with self.assertRaisesRegex(tarfile.ReadError, 'Tarfile is invalid.'):
            fetch_card_data(version='3.2.6')

    def test_extract_sequence_dna(self):
        with open(self.get_data_path('card_test_data/card_test.json'), 'rb') as f:
            db = json.load(f)
        obs = extract_sequence('dna', '2', '1188', db)
        exp = DNA("ATGAAAGCATATTTCATCGCCATACTTACCTTATTCACTTGTATAGCTACCGTCGTCCGGGCGCAGCAAATGTCTGAACTTGAAAACCGGATTGACAGTCTGCTCAATGGCAAGAAAGCCACCGTTGGTATAGCCGTATGGACAGACAAAGGAGACATGCTCCGGTATAACGACCATGTACACTTCCCCTTGCTCAGTGTATTCAAATTCCATGTGGCACTGGCCGTACTGGACAAGATGGATAAGCAAAGCATCAGTCTGGACAGCATTGTTTCCATAAAGGCATCCCAAATGCCGCCCAATACCTACAGCCCCCTGCGGAAGAAGTTTCCCGACCAGGATTTCACGATTACGCTTAGGGAACTGATGCAATACAGCATTTCCCAAAGCGACAACAATGCCTGCGACATCTTGATAGAATATGCAGGAGGCATCAAACATATCAACGACTATATCCACCGGTTGAGTATCGACTCCTTCAACCTCTCGGAAACAGAAGACGGCATGCACTCCAGCTTCGAGGCTGTATACCGCAACTGGAGTACTCCTTCCGCTATGGTCCGACTACTGAGAACGGCTGATGAAAAAGAGTTGTTCTCCAACAAGGAGCTGAAAGACTTCTTGTGGCAGACCATGATAGATACTGAAACCGGTGCCAACAAACTGAAAGGTATGTTGCCAGCCAAAACCGTGGTAGGACACAAGACCGGCTCTTCCGACCGCAATGCCGACGGTATGAAAACTGCAGATAATGATGCCGGCCTCGTTATCCTTCCCGACGGCCGGAAATACTACATTGCCGCCTTCGTCATGGACTCATACGAGACGGATGAGGACAATGCGAACATCATCGCCCGCATATCACGCATGGTATATGATGCGATGAGATGA")
        exp.metadata['id'] = 'gb|GQ343019.1|+|132-1023|ARO:3002999|CblA-1'
        exp.metadata['description'] = '[mixed culture bacterium AX_gF3SD01_15]'
        self.assertEqual(exp, obs)
        self.assertIsInstance(obs, DNA)

    def test_extract_sequence_protein(self):
        with open(self.get_data_path('card_test_data/card_test.json'), 'rb') as f:
            db = json.load(f)
        obs = extract_sequence('protein', '2', '1188', db)
        exp = Protein("MKAYFIAILTLFTCIATVVRAQQMSELENRIDSLLNGKKATVGIAVWTDKGDMLRYNDHVHFPLLSVFKFHVALAVLDKMDKQSISLDSIVSIKASQMPPNTYSPLRKKFPDQDFTITLRELMQYSISQSDNNACDILIEYAGGIKHINDYIHRLSIDSFNLSETEDGMHSSFEAVYRNWSTPSAMVRLLRTADEKELFSNKELKDFLWQTMIDTETGANKLKGMLPAKTVVGHKTGSSDRNADGMKTADNDAGLVILPDGRKYYIAAFVMDSYETDEDNANIIARISRMVYDAMR")
        exp.metadata['id'] = 'gb|ACT97415.1|ARO:3002999|CblA-1'
        exp.metadata['description'] = '[mixed culture bacterium AX_gF3SD01_15]'
        self.assertEqual(exp, obs)
        self.assertIsInstance(obs, Protein)

    def test_read_from_card_file(self):
        path = self.get_data_path('card.json')
        genomes = _read_from_card_file(path, 'protein')
        self.assertIsInstance(genomes, ProteinFASTAFormat)


    def test_read_from_card_generator(self):
        path = self.get_data_path('card.json')
        genomes = _read_from_card_file(path, 'protein')
        generator = ProteinIterator(genomes)
        self.assertIsInstance(generator, ProteinIterator)


class TestCARDAnnotationTypesAndFormats(AMRTypesTestPluginBase):

    def test_df_to_card_annotation_format_transformer(self):
        filepath = self.get_data_path('rgi_output.txt')
        transformer = self.get_transformer(pd.DataFrame, CARDAnnotationFormat)
        df = pd.read_csv(filepath, sep="\t")
        obs = transformer(df)
        self.assertIsInstance(obs, CARDAnnotationFormat)

    def test_card_annotation_format_to_df_transformer(self):
        filepath = self.get_data_path('rgi_output.txt')
        transformer = self.get_transformer(CARDAnnotationFormat, pd.DataFrame)
        card_anno = CARDAnnotationFormat(filepath, mode='r')
        obs = transformer(card_anno)
        self.assertIsInstance(obs, pd.DataFrame)

    def test_card_annotation(self):
        exp = pd.read_csv("/Users/vinzent/Desktop/bokulich_project/Code/RESCRIPt/rescript/types/tests/data/rgi_output.txt", sep='\t')
        obs = card_annotation(input_seq='/Users/vinzent/Desktop/bokulich_project/Code/RESCRIPt/rescript/types/tests/data/rgi_input.fna', output='test')
        assert_frame_equal(exp, obs)

    def test_card_annotation_txt_to_fasta(self):
        filepath = self.get_data_path('rgi_output.txt')
        df = pd.read_csv(filepath, sep="\t")
        protein = card_annotation_txt_to_fasta(df)[0]
        dna = card_annotation_txt_to_fasta(df)[1]
        format = ProteinFASTAFormat(str(protein), mode='r')
        print(format)
        print(dna)
