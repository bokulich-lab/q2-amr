import os
import shutil
import subprocess
import tempfile
from unittest.mock import patch

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from qiime2.plugin.testing import TestPluginBase

from q2_amr.card import (
    annotate_reads_card,
    create_count_table,
    extract_sample_stats,
    move_files,
    plot_sample_stats,
    read_in_txt,
    run_rgi_bwt,
    visualize_annotation_stats,
)
from q2_amr.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDDatabaseFormat,
    CARDGeneAnnotationDirectoryFormat,
)


class TestAnnotateReadsCARD(TestPluginBase):
    package = "q2_amr.tests"

    def test_annotate_reads_card_single(self):
        self.annotate_reads_card_test_body("single")

    def test_annotate_reads_card_paired(self):
        self.annotate_reads_card_test_body("paired")

    def annotate_reads_card_test_body(self, read_type):
        manifest = self.get_data_path(f"MANIFEST_reads_{read_type}")
        if read_type == "single":
            reads = SingleLanePerSampleSingleEndFastqDirFmt()
        else:
            reads = SingleLanePerSamplePairedEndFastqDirFmt()
        shutil.copy(manifest, os.path.join(str(reads), "MANIFEST"))
        output_allele = self.get_data_path("allele_mapping_data.txt")
        output_gene = self.get_data_path("gene_mapping_data.txt")
        output_stats = self.get_data_path("overall_mapping_stats.txt")
        card_db = CARDDatabaseFormat()

        def mock_run_rgi_bwt(cwd, samp, **kwargs):
            shutil.copy(output_allele, f"{cwd}/{samp}/allele_mapping_data.txt")
            shutil.copy(output_gene, f"{cwd}/{samp}/gene_mapping_data.txt")
            shutil.copy(output_stats, f"{cwd}/{samp}/overall_mapping_stats.txt")

        with patch("q2_amr.card.run_rgi_bwt", side_effect=mock_run_rgi_bwt), patch(
            "q2_amr.card.load_preprocess_card_db"
        ):
            result = annotate_reads_card(reads, card_db)
            self.assertIsInstance(result[0], CARDAlleleAnnotationDirectoryFormat)
            self.assertIsInstance(result[1], CARDGeneAnnotationDirectoryFormat)
            self.assertIsInstance(result[2], pd.DataFrame)
            self.assertIsInstance(result[3], pd.DataFrame)
            self.assertTrue(
                os.path.exists(
                    os.path.join(str(result[0]), "sample1", "allele_mapping_data.txt")
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(str(result[1]), "sample1", "gene_mapping_data.txt")
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(str(result[0]), "sample1", "overall_mapping_stats.txt")
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(str(result[1]), "sample1", "overall_mapping_stats.txt")
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(str(result[0]), "sample2", "allele_mapping_data.txt")
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(str(result[1]), "sample2", "gene_mapping_data.txt")
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(str(result[0]), "sample2", "overall_mapping_stats.txt")
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(str(result[1]), "sample1", "overall_mapping_stats.txt")
                )
            )

    def test_run_rgi_bwt(self):
        with patch("q2_amr.card.run_command") as mock_run_command:
            run_rgi_bwt(
                "path_tmp",
                "sample1",
                "path_fwd",
                "path_rev",
                "bowtie2",
                8,
                True,
                3,
                5,
                3.2,
            )
            mock_run_command.assert_called_once_with(
                [
                    "rgi",
                    "bwt",
                    "--read_one",
                    "path_fwd",
                    "--output_file",
                    "path_tmp/sample1/sample1",
                    "-n",
                    "8",
                    "--local",
                    "--clean",
                    "--aligner",
                    "bowtie2",
                    "--read_two",
                    "path_rev",
                    "--include_baits",
                    "--mapq",
                    "3",
                    "--mapped",
                    "5",
                    "--coverage",
                    "3.2",
                ],
                "path_tmp",
                verbose=True,
            )

    @patch("q2_amr.card.run_command")
    def test_exception_raised(self, mock_run_command):
        mock_run_command.side_effect = subprocess.CalledProcessError(1, "cmd")
        cwd = "path/cwd"
        samp = "sample1"
        fwd = "path/fwd"
        rev = "path/rev"
        aligner = "bwa"
        threads = 1
        include_baits = True
        mapq = 0.3
        mapped = 0.3
        coverage = 0.3
        expected_message = (
            "An error was encountered while running rgi, "
            "(return code 1), please inspect stdout and stderr to learn more."
        )
        with self.assertRaises(Exception) as cm:
            run_rgi_bwt(
                cwd,
                samp,
                fwd,
                rev,
                aligner,
                threads,
                include_baits,
                mapq,
                mapped,
                coverage,
            )
        self.assertEqual(str(cm.exception), expected_message)

    def test_move_files_allele(self):
        self.move_files_test_body("allele")

    def test_move_files_gene(self):
        self.move_files_test_body("gene")

    def move_files_test_body(self, map_type):
        with tempfile.TemporaryDirectory() as tmp:
            samp = "sample1"
            des_dir = os.path.join(tmp, "des_dir")
            os.makedirs(os.path.join(tmp, samp))
            os.makedirs(os.path.join(des_dir, samp))
            with open(
                os.path.join(tmp, samp, f"{map_type}_mapping_data.txt"), "w"
            ) as file:
                file.write("Sample mapping data")
            with open(
                os.path.join(tmp, samp, "overall_mapping_stats.txt"), "w"
            ) as file:
                file.write("Overall mapping stats")
            move_files(tmp, des_dir, samp, map_type)
            self.assertFalse(
                os.path.exists(os.path.join(tmp, samp, f"{map_type}_mapping_data.txt"))
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        des_dir,
                        samp,
                        f"{map_type}_mapping_data.txt",
                    )
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        des_dir,
                        samp,
                        "overall_mapping_stats.txt",
                    )
                )
            )

    def test_extract_sample_stats(self):
        with tempfile.TemporaryDirectory() as tmp:
            samp = "sample1"
            sample_stats = {}
            mapping_stats_path = self.get_data_path("overall_mapping_stats.txt")
            shutil.copy(mapping_stats_path, tmp)
            extract_sample_stats(tmp, samp, sample_stats)
            expected_result = {
                "sample1": {
                    "total_reads": 5000,
                    "mapped_reads": 59,
                    "percentage": 1.18,
                }
            }
            self.assertEqual(sample_stats, expected_result)

    def test_plot_sample_stats(self):
        with tempfile.TemporaryDirectory() as tmp:
            sample_stats = {
                "sample1": {
                    "total_reads": 5000,
                    "mapped_reads": 106,
                    "percentage": 2.12,
                },
                "sample2": {
                    "total_reads": 7000,
                    "mapped_reads": 212,
                    "percentage": 3.03,
                },
            }
            plot_sample_stats(sample_stats, tmp)
            self.assertTrue(os.path.exists(os.path.join(tmp, "sample_stats_plot.html")))

    def test_visualize_annotation_stats(self):
        def mock_extract_sample_stats(samp_dir, samp, sample_stats):
            sample_stats[samp] = {
                "total_reads": 5000,
                "mapped_reads": 106,
                "percentage": 2.12,
            }

        def mock_plot_sample_stats(sample_stats, output_dir):
            with open(os.path.join(tmp, "sample_stats_plot.html"), "w") as file:
                file.write("file")

        amr_reads_annotation = CARDGeneAnnotationDirectoryFormat()
        sample1_dir = os.path.join(str(amr_reads_annotation), "sample1")
        sample2_dir = os.path.join(str(amr_reads_annotation), "sample2")
        os.makedirs(sample1_dir)
        os.makedirs(sample2_dir)
        with patch(
            "q2_amr.card.extract_sample_stats", side_effect=mock_extract_sample_stats
        ), patch(
            "q2_amr.card.plot_sample_stats", side_effect=mock_plot_sample_stats
        ), tempfile.TemporaryDirectory() as tmp:
            visualize_annotation_stats(tmp, amr_reads_annotation)
            self.assertTrue(os.path.exists(os.path.join(tmp, "sample_stats_plot.html")))
            self.assertTrue(os.path.exists(os.path.join(tmp, "index.html")))
            self.assertTrue(os.path.exists(os.path.join(tmp, "q2templateassets")))

    mapping_data_sample1 = pd.DataFrame(
        {
            "ARO Accession": [3000796, 3000815, 3000805, 3000026],
            "sample1": [1, 1, 1, 1],
        }
    )

    mapping_data_sample2 = pd.DataFrame(
        {
            "ARO Accession": [3000797, 3000815, 3000805, 3000026],
            "sample2": [1, 1, 1, 2],
        }
    )

    count_table = pd.DataFrame(
        {
            "index": ["sample1", "sample2"],
            3000796: [1, 0],
            3000815: [1, 1],
            3000805: [1, 1],
            3000026: [1, 2],
            3000797: [0, 1],
        }
    )

    def test_read_in_txt_allele(self):
        self.read_in_txt_test_body("allele", self.mapping_data_sample1)

    def test_read_in_txt_gene(self):
        self.read_in_txt_test_body("gene", self.mapping_data_sample1)

    def read_in_txt_test_body(self, map_type, mapping_data):
        samp = "sample1"
        mapping_file = self.get_data_path(f"{map_type}_mapping_data.txt")
        exp = mapping_data
        with tempfile.TemporaryDirectory() as tmp:
            samp_dir = os.path.join(tmp, samp)
            os.mkdir(samp_dir)
            shutil.copy(mapping_file, samp_dir)
            obs = read_in_txt(tmp, samp, map_type)
            obs["ARO Accession"] = obs["ARO Accession"].astype(int)
            pd.testing.assert_frame_equal(exp, obs)

    def test_create_count_table(self):
        df_list = [self.mapping_data_sample1, self.mapping_data_sample2]
        obs = create_count_table(df_list)
        exp = self.count_table
        exp.set_index("index", inplace=True)
        exp.index.name = "sample_id"
        exp = exp.astype(float)
        exp.columns = exp.columns.astype(float)
        pd.testing.assert_frame_equal(exp, obs)
