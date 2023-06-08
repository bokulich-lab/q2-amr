import os
import shutil
import subprocess
import tempfile
from copy import deepcopy
from unittest.mock import ANY, MagicMock, call, patch

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
            shutil.copy(manifest, os.path.join(str(reads), "MANIFEST"))

        else:
            reads = SingleLanePerSamplePairedEndFastqDirFmt()
            shutil.copy(manifest, os.path.join(str(reads), "MANIFEST"))

        output_allele = self.get_data_path("output.allele_mapping_data.txt")
        output_gene = self.get_data_path("output.gene_mapping_data.txt")
        output_stats = self.get_data_path("output.overall_mapping_stats.txt")
        card_db = CARDDatabaseFormat()

        def copy_needed_files(cwd, samp, **kwargs):
            samp_dir = os.path.join(cwd, samp)
            shutil.copy(output_allele, samp_dir)
            shutil.copy(output_gene, samp_dir)
            shutil.copy(output_stats, samp_dir)

        def return_count_table(df_list):
            count_table = deepcopy(self.table)
            count_table.set_index("sample_id", inplace=True)
            count_table = count_table.astype(float)
            count_table.columns = count_table.columns.astype(float)
            return count_table

        mock_run_rgi_bwt = MagicMock(side_effect=copy_needed_files)
        mock_run_rgi_load = MagicMock()
        mock_read_in_txt = MagicMock()
        mock_create_count_table = MagicMock(side_effect=return_count_table)
        with patch("q2_amr.card.run_rgi_bwt", mock_run_rgi_bwt), patch(
            "q2_amr.card.load_preprocess_card_db", mock_run_rgi_load
        ), patch("q2_amr.card.read_in_txt", mock_read_in_txt), patch(
            "q2_amr.card.create_count_table", mock_create_count_table
        ):
            result = annotate_reads_card(reads, card_db)
            first_call_args = mock_run_rgi_bwt.call_args_list[0]
            tmp_dir = first_call_args.kwargs["cwd"]
            if read_type == "single":
                exp_calls_mock_run = [
                    call(
                        cwd=tmp_dir,
                        samp="sample1",
                        fwd=f"{reads}/sample1_00_L001_R1_001.fastq.gz",
                        aligner="kma",
                        rev=None,
                        threads=1,
                        include_baits=False,
                        mapq=None,
                        mapped=None,
                        coverage=None,
                    ),
                    call(
                        cwd=tmp_dir,
                        samp="sample2",
                        fwd=f"{reads}/sample2_00_L001_R1_001.fastq.gz",
                        aligner="kma",
                        rev=None,
                        threads=1,
                        include_baits=False,
                        mapq=None,
                        mapped=None,
                        coverage=None,
                    ),
                ]
            else:
                exp_calls_mock_run = [
                    call(
                        cwd=tmp_dir,
                        samp="sample1",
                        fwd=f"{reads}/sample1_00_L001_R1_001.fastq.gz",
                        rev=f"{reads}/sample1_00_L001_R2_001.fastq.gz",
                        aligner="kma",
                        threads=1,
                        include_baits=False,
                        mapq=None,
                        mapped=None,
                        coverage=None,
                    ),
                    call(
                        cwd=tmp_dir,
                        samp="sample2",
                        fwd=f"{reads}/sample2_00_L001_R1_001.fastq.gz",
                        rev=f"{reads}/sample2_00_L001_R2_001.fastq.gz",
                        aligner="kma",
                        threads=1,
                        include_baits=False,
                        mapq=None,
                        mapped=None,
                        coverage=None,
                    ),
                ]
            exp_calls_mock_load = [
                call(tmp_dir, ANY, "load"),
                call(tmp_dir, ANY, "preprocess"),
                call(tmp_dir, ANY, "load_fasta"),
            ]
            exp_calls_mock_read = [
                call(f"{tmp_dir}/sample1", "allele"),
                call(f"{tmp_dir}/sample1", "gene"),
                call(f"{tmp_dir}/sample2", "allele"),
                call(f"{tmp_dir}/sample2", "gene"),
            ]
            exp_calls_mock_count = [call([ANY, ANY]), call([ANY, ANY])]
            mock_run_rgi_bwt.assert_has_calls(exp_calls_mock_run)
            mock_run_rgi_load.assert_has_calls(exp_calls_mock_load)
            mock_read_in_txt.assert_has_calls(exp_calls_mock_read)
            mock_create_count_table.assert_has_calls(exp_calls_mock_count)
            self.assertIsInstance(result[0], CARDAlleleAnnotationDirectoryFormat)
            self.assertIsInstance(result[1], CARDGeneAnnotationDirectoryFormat)
            self.assertIsInstance(result[2], pd.DataFrame)
            self.assertIsInstance(result[3], pd.DataFrame)
            for num in [0, 1]:
                map_type = "allele" if num == 0 else "gene"
                for samp in ["sample1", "sample2"]:
                    for file in [
                        f"{map_type}_mapping_data.txt",
                        "overall_mapping_stats.txt",
                    ]:
                        self.assertTrue(
                            os.path.exists(os.path.join(str(result[num]), samp, file))
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
                    "path_tmp/sample1/output",
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
        expected_message = (
            "An error was encountered while running rgi, "
            "(return code 1), please inspect stdout and stderr to learn more."
        )
        with self.assertRaises(Exception) as cm:
            run_rgi_bwt(
                cwd="path/cwd",
                samp="sample1",
                fwd="path/fwd",
                rev="path/rev",
                aligner="bwa",
                threads=1,
                include_baits=True,
                mapq=0.3,
                mapped=0.3,
                coverage=0.3,
            )
        self.assertEqual(str(cm.exception), expected_message)

    def test_move_files_allele(self):
        self.move_files_test_body("allele")

    def test_move_files_gene(self):
        self.move_files_test_body("gene")

    def move_files_test_body(self, map_type):
        with tempfile.TemporaryDirectory() as tmp:
            source_dir = os.path.join(tmp, "source_dir")
            des_dir = os.path.join(
                tmp,
                "des_dir",
            )
            os.makedirs(os.path.join(source_dir))
            os.makedirs(os.path.join(des_dir))
            mapping_data = self.get_data_path(f"output.{map_type}_mapping_data.txt")
            mapping_stats = self.get_data_path("output.overall_mapping_stats.txt")
            shutil.copy(mapping_data, source_dir)
            shutil.copy(mapping_stats, source_dir)
            move_files(source_dir, des_dir, map_type)
            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        des_dir,
                        f"{map_type}_mapping_data.txt",
                    )
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        des_dir,
                        "overall_mapping_stats.txt",
                    )
                )
            )

    def test_extract_sample_stats(self):
        with tempfile.TemporaryDirectory() as tmp:
            mapping_stats_path = self.get_data_path("output.overall_mapping_stats.txt")
            new_mapping_stats_path = os.path.join(tmp, "overall_mapping_stats.txt")
            shutil.copy(mapping_stats_path, new_mapping_stats_path)
            sample_stats = extract_sample_stats(tmp)
            expected_result = {
                "total_reads": 5000,
                "mapped_reads": 59,
                "percentage": 1.18,
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
        def mock_plot_sample_stats(sample_stats, output_dir):
            with open(os.path.join(tmp, "sample_stats_plot.html"), "w") as file:
                file.write("file")

        amr_reads_annotation = CARDGeneAnnotationDirectoryFormat()
        sample1_dir = os.path.join(str(amr_reads_annotation), "sample1")
        sample2_dir = os.path.join(str(amr_reads_annotation), "sample2")
        os.makedirs(sample1_dir)
        os.makedirs(sample2_dir)
        with patch("q2_amr.card.extract_sample_stats"), patch(
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

    table = pd.DataFrame(
        {
            "sample_id": ["sample1", "sample2"],
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
        mapping_file = self.get_data_path(f"output.{map_type}_mapping_data.txt")
        exp = mapping_data
        with tempfile.TemporaryDirectory() as tmp:
            samp_dir = os.path.join(tmp, "sample1")
            os.mkdir(samp_dir)
            shutil.copy(mapping_file, samp_dir)
            obs = read_in_txt(samp_dir, map_type)
            obs["ARO Accession"] = obs["ARO Accession"].astype(int)
            pd.testing.assert_frame_equal(exp, obs)

    def test_create_count_table(self):
        df_list = [self.mapping_data_sample1, self.mapping_data_sample2]
        obs = create_count_table(df_list)
        exp = self.table
        exp.set_index("sample_id", inplace=True)
        exp = exp.astype(float)
        exp.columns = exp.columns.astype(float)
        pd.testing.assert_frame_equal(exp, obs)
