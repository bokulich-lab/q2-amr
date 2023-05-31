import os
import shutil
import tempfile
from unittest.mock import patch

import pandas as pd
from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt, \
    SingleLanePerSampleSingleEndFastqDirFmt
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase

from q2_amr.card import load_card_db, preprocess_card_db, load_card_db_fasta, annotate_reads_card
from q2_amr.types import CARDDatabaseFormat, CARDAlleleAnnotationDirectoryFormat, CARDGeneAnnotationDirectoryFormat


class TestBwt(TestPluginBase):
    package = 'q2_amr.tests'

    def test_load_card_db(self):
        card_db = Artifact.load(self.get_data_path('card_db_test.qza'))
        a = card_db.view(CARDDatabaseFormat)
        with tempfile.TemporaryDirectory() as tmp:
            load_card_db(tmp, a)
            assert os.path.exists(os.path.join(tmp, 'localDB', 'card.json'))
            assert os.path.exists(os.path.join(tmp, 'localDB', 'loaded_databases.json'))
            preprocess_card_db(tmp, a)
            assert os.path.exists(os.path.join(tmp, 'card_database_v3.2.5.fasta'))
            load_card_db_fasta(tmp, a)
            assert os.path.exists(os.path.join(tmp, 'localDB', 'card_reference.fasta'))

    def test_preprocess_card_db(self):
        card_db = Artifact.load(self.get_data_path('card_db_test.qza'))
        a = card_db.view(CARDDatabaseFormat)
        with tempfile.TemporaryDirectory() as tmp:
            preprocess_card_db(tmp, a)
            assert os.path.exists(os.path.join(tmp, 'card_database_v3.2.5.fasta'))

    def test_annotate_reads_card(self):
        output_allele = self.get_data_path('sample1.allele_mapping_data.txt')
        output_gene = self.get_data_path('sample1.gene_mapping_data.txt')
        output_stats = self.get_data_path('sample1.overall_mapping_stats.txt')
        reads = SingleLanePerSampleSingleEndFastqDirFmt()
        card_db = CARDDatabaseFormat()
        manifest = self.get_data_path('MANIFEST_reads')
        shutil.copy(manifest, os.path.join(str(reads), 'MANIFEST'))

        def mock_run_rgi_bwt(tmp, samp, fwd, rev, aligner, threads, include_baits, mapq, mapped, coverage):
            shutil.copy(output_allele, f"{tmp}/{samp}/{samp}.allele_mapping_data.txt")
            shutil.copy(output_gene, f"{tmp}/{samp}/{samp}.gene_mapping_data.txt")
            shutil.copy(output_stats, f"{tmp}/{samp}/{samp}.overall_mapping_stats.txt")

        with patch('q2_amr.card.run_rgi_bwt', side_effect=mock_run_rgi_bwt), patch('q2_amr.card.load_card_db'), \
             patch('q2_amr.card.preprocess_card_db'), patch('q2_amr.card.load_card_db_fasta'):
            result = annotate_reads_card(reads, card_db)
            print('a')
            self.assertIsInstance(result[0], CARDAlleleAnnotationDirectoryFormat)
            self.assertIsInstance(result[1], CARDGeneAnnotationDirectoryFormat)
            self.assertIsInstance(result[2], pd.DataFrame)
            self.assertIsInstance(result[3], pd.DataFrame)
            self.assertTrue(os.path.exists(os.path.join(str(result[0]), 'sample1', 'sample1.allele_mapping_data.txt')))
            self.assertTrue(os.path.exists(os.path.join(str(result[1]), 'sample1', 'sample1.gene_mapping_data.txt')))
            self.assertTrue(
                os.path.exists(os.path.join(str(result[0]), 'sample1', 'sample1.overall_mapping_stats.txt')))
            self.assertTrue(
                os.path.exists(os.path.join(str(result[1]), 'sample1', 'sample1.overall_mapping_stats.txt')))
            self.assertTrue(os.path.exists(os.path.join(str(result[0]), 'sample2', 'sample2.allele_mapping_data.txt')))
            self.assertTrue(os.path.exists(os.path.join(str(result[1]), 'sample2', 'sample2.gene_mapping_data.txt')))
            self.assertTrue(
                os.path.exists(os.path.join(str(result[0]), 'sample2', 'sample2.overall_mapping_stats.txt')))
            self.assertTrue(
                os.path.exists(os.path.join(str(result[1]), 'sample1', 'sample1.overall_mapping_stats.txt')))
