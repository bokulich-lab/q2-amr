import os
import shutil
import tempfile
from unittest.mock import patch
import pandas as pd
from q2_types.per_sample_sequences import SingleLanePerSamplePairedEndFastqDirFmt, \
    SingleLanePerSampleSingleEndFastqDirFmt
from qiime2.plugin.testing import TestPluginBase
from q2_amr.card import load_card_db, preprocess_card_db, load_card_db_fasta, annotate_reads_card, run_rgi_bwt, \
    move_files, extract_sample_stats, plot_sample_stats, visualize_annotation_stats
from q2_amr.types import CARDDatabaseFormat, CARDAlleleAnnotationDirectoryFormat, CARDGeneAnnotationDirectoryFormat


class TestAnnotateReadsCARD(TestPluginBase):
    package = 'q2_amr.tests'

    def test_load_card_db(self):
        card_db = CARDDatabaseFormat()
        with patch('q2_amr.card.run_command') as mock_run_command:
            load_card_db('path_tmp', card_db)
            mock_run_command.assert_called_once_with(['rgi', 'load', '--card_json', str(card_db), '--local'],
                                                     'path_tmp', verbose=True)

    def test_preprocess_card_db(self):
        card_db = CARDDatabaseFormat()
        with patch('q2_amr.card.run_command') as mock_run_command:
            preprocess_card_db('path_tmp', card_db)
            mock_run_command.assert_called_once_with(['rgi', 'card_annotation', '-i', str(card_db)],
                                                     'path_tmp', verbose=True)

    def test_load_card_db_fasta(self):
        card_db = self.get_data_path('card_test.json')
        with patch('q2_amr.card.run_command') as mock_run_command:
            load_card_db_fasta('path_tmp', card_db)
            mock_run_command.assert_called_once_with(['rgi', 'load', '-i', str(card_db), '--card_annotation',
                                                      f'card_database_v3.2.5.fasta', '--local'], 'path_tmp',
                                                     verbose=True)

    def test_annotate_reads_card_single(self):
        manifest = self.get_data_path('MANIFEST_reads_single')
        reads = SingleLanePerSampleSingleEndFastqDirFmt()
        shutil.copy(manifest, os.path.join(str(reads), 'MANIFEST'))
        self.annotate_reads_card_test_body(reads)

    def test_annotate_reads_card_paired(self):
        manifest = self.get_data_path('MANIFEST_reads_paired')
        reads = SingleLanePerSamplePairedEndFastqDirFmt()
        shutil.copy(manifest, os.path.join(str(reads), 'MANIFEST'))
        self.annotate_reads_card_test_body(reads)

    def annotate_reads_card_test_body(self, reads):
        output_allele = self.get_data_path('sample1.allele_mapping_data.txt')
        output_gene = self.get_data_path('sample1.gene_mapping_data.txt')
        output_stats = self.get_data_path('sample1.overall_mapping_stats.txt')
        card_db = CARDDatabaseFormat()

        def mock_run_rgi_bwt(tmp, samp, fwd, rev, aligner, threads, include_baits, mapq, mapped, coverage):
            shutil.copy(output_allele, f"{tmp}/{samp}/{samp}.allele_mapping_data.txt")
            shutil.copy(output_gene, f"{tmp}/{samp}/{samp}.gene_mapping_data.txt")
            shutil.copy(output_stats, f"{tmp}/{samp}/{samp}.overall_mapping_stats.txt")

        with patch('q2_amr.card.run_rgi_bwt', side_effect=mock_run_rgi_bwt), patch('q2_amr.card.load_card_db'), \
                patch('q2_amr.card.preprocess_card_db'), patch('q2_amr.card.load_card_db_fasta'):
            result = annotate_reads_card(reads, card_db)
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

    def test_run_rgi_bwt(self):
        with patch('q2_amr.card.run_command') as mock_run_command:
            run_rgi_bwt('path_tmp', 'sample1', 'path_fwd', 'path_rev', 'bowtie2', 8, True, 3, 5, 3.2)
            mock_run_command.assert_called_once_with(['rgi', 'bwt', '--read_one', 'path_fwd', '--output_file',
                                                      'path_tmp/sample1/sample1', '-n', '8', '--local', '--clean',
                                                      '--read_two', 'path_rev', '--include_baits', '--mapq', '3',
                                                      '--mapped', '5', '--coverage', '3.2', '--aligner', 'bowtie2'],
                                                     'path_tmp', verbose=True)

    def test_move_files_allele(self):
        map_type = "allele"
        des_dir = CARDAlleleAnnotationDirectoryFormat()
        self.move_files_test_body(map_type, des_dir)

    def test_move_files_gene(self):
        map_type = "gene"
        des_dir = CARDGeneAnnotationDirectoryFormat()
        self.move_files_test_body(map_type, des_dir)

    def move_files_test_body(self, map_type, des_dir):
        with tempfile.TemporaryDirectory() as tmp:
            samp = "sample1"
            os.makedirs(os.path.join(tmp, samp))
            os.makedirs(os.path.join(str(des_dir), samp))
            with open(os.path.join(tmp, samp, f"{samp}.{map_type}_mapping_data.txt"), "w") as file:
                file.write("Sample mapping data")
            with open(os.path.join(tmp, samp, f"{samp}.overall_mapping_stats.txt"), "w") as file:
                file.write("Overall mapping stats")
            move_files(tmp, samp, map_type, des_dir)
            self.assertFalse(os.path.exists(os.path.join(tmp, samp, f"{samp}.{map_type}_mapping_data.txt")))
            self.assertTrue(
                os.path.exists(os.path.join(os.path.join(str(des_dir), samp), f"{samp}.{map_type}_mapping_data.txt")))
            self.assertTrue(
                os.path.exists(os.path.join(os.path.join(str(des_dir), samp), f"{samp}.overall_mapping_stats.txt")))

    def test_extract_sample_stats(self):
        with tempfile.TemporaryDirectory() as tmp:
            samp = "sample1"
            sample_stats = {}
            stats_content = """
            **********************************************
            Stats for BAM file(s): 
            **********************************************
    
            Total reads:       5000
            Mapped reads:      106	(2.12%)
            Forward strand:    4947	(98.94%)
            Reverse strand:    53	(1.06%)
            Failed QC:         0	(0%)
            Duplicates:        0	(0%)
            Paired-end reads:  5000	(100%)
            'Proper-pairs':    22	(0.44%)
            Both pairs mapped: 76	(1.52%)
            Read 1:            2500
            Read 2:            2500
            Singletons:        30	(0.6%)
            """
            with open(os.path.join(tmp, f"{samp}.overall_mapping_stats.txt"), "w") as file:
                file.write(stats_content)
            extract_sample_stats(tmp, samp, sample_stats)
            expected_result = {'sample1': {'total_reads': 5000, 'mapped_reads': 106, 'percentage': 2.12}}
            self.assertEqual(sample_stats, expected_result)

    def test_plot_sample_stats(self):
        with tempfile.TemporaryDirectory() as tmp:
            sample_stats = {
                'sample1': {'total_reads': 5000, 'mapped_reads': 106, 'percentage': 2.12},
                'sample2': {'total_reads': 7000, 'mapped_reads': 212, 'percentage': 3.03}
            }
            plot_sample_stats(sample_stats, tmp)
            self.assertTrue(os.path.exists(os.path.join(tmp, "sample_stats_plot.html")))

    def test_visualize_annotation_stats(self):
        def mock_extract_sample_stats(samp_dir, samp, sample_stats):
            sample_stats[samp] = {'total_reads': 5000, 'mapped_reads': 106, 'percentage': 2.12}

        def mock_plot_sample_stats(sample_stats, output_dir):
            with open(os.path.join(tmp, "sample_stats_plot.html"), "w") as file:
                file.write("file")

        sample_stats = {}
        amr_reads_annotation = CARDGeneAnnotationDirectoryFormat()
        sample1_dir = os.path.join(str(amr_reads_annotation), "sample1")
        sample2_dir = os.path.join(str(amr_reads_annotation), "sample2")
        os.makedirs(sample1_dir)
        os.makedirs(sample2_dir)
        with patch('q2_amr.card.extract_sample_stats', side_effect=mock_extract_sample_stats), \
                patch('q2_amr.card.plot_sample_stats', side_effect=mock_plot_sample_stats), \
                tempfile.TemporaryDirectory() as tmp:
            visualize_annotation_stats(tmp, amr_reads_annotation)
            self.assertTrue(os.path.exists(os.path.join(tmp, "sample_stats_plot.html")))
            self.assertTrue(os.path.exists(os.path.join(tmp, "index.html")))
            self.assertTrue(os.path.exists(os.path.join(tmp, "q2templateassets")))
