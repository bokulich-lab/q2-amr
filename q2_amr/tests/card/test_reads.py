import os
import shutil
import subprocess
import tempfile
from unittest.mock import ANY, MagicMock, call, patch

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from qiime2.plugin.testing import TestPluginBase
from test_mags import TestAnnotateMagsCard

from q2_amr.card.reads import (
    annotate_reads_card,
    extract_sample_stats,
    move_files,
    plot_sample_stats,
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

    def copy_needed_files(self, cwd, samp, **kwargs):
        output_allele = self.get_data_path("output.allele_mapping_data.txt")
        output_gene = self.get_data_path("output.gene_mapping_data.txt")
        output_stats = self.get_data_path("output.overall_mapping_stats.txt")
        samp_dir = os.path.join(cwd, samp)
        shutil.copy(output_allele, samp_dir)
        shutil.copy(output_gene, samp_dir)
        shutil.copy(output_stats, samp_dir)

    def annotate_reads_card_test_body(self, read_type):
        manifest = self.get_data_path(f"MANIFEST_reads_{read_type}")
        if read_type == "single":
            reads = SingleLanePerSampleSingleEndFastqDirFmt()
            shutil.copy(manifest, os.path.join(str(reads), "MANIFEST"))
        else:
            reads = SingleLanePerSamplePairedEndFastqDirFmt()
            shutil.copy(manifest, os.path.join(str(reads), "MANIFEST"))
        card_db = CARDDatabaseFormat()
        mock_run_rgi_bwt = MagicMock(side_effect=self.copy_needed_files)
        mock_run_rgi_load = MagicMock()
        mock_read_in_txt = MagicMock()
        mag_test_class = TestAnnotateMagsCard()
        mock_create_count_table = MagicMock(
            side_effect=mag_test_class.return_count_table
        )
        with patch("q2_amr.card.reads.run_rgi_bwt", mock_run_rgi_bwt), patch(
            "q2_amr.card.reads.load_card_db", mock_run_rgi_load
        ), patch("q2_amr.card.reads.read_in_txt", mock_read_in_txt), patch(
            "q2_amr.card.reads.create_count_table", mock_create_count_table
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
                        include_wildcard=False,
                        include_other_models=False,
                    ),
                    call(
                        cwd=tmp_dir,
                        samp="sample2",
                        fwd=f"{reads}/sample2_00_L001_R1_001.fastq.gz",
                        aligner="kma",
                        rev=None,
                        threads=1,
                        include_wildcard=False,
                        include_other_models=False,
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
                        include_wildcard=False,
                        include_other_models=False,
                    ),
                    call(
                        cwd=tmp_dir,
                        samp="sample2",
                        fwd=f"{reads}/sample2_00_L001_R1_001.fastq.gz",
                        rev=f"{reads}/sample2_00_L001_R2_001.fastq.gz",
                        aligner="kma",
                        threads=1,
                        include_wildcard=False,
                        include_other_models=False,
                    ),
                ]
            exp_calls_mock_load = [
                call(
                    tmp=tmp_dir,
                    card_db=ANY,
                    fasta=True,
                    include_other_models=False,
                    include_wildcard=False,
                ),
            ]
            exp_calls_mock_read = [
                call(
                    path=f"{tmp_dir}/sample1/output.allele_mapping_data.txt",
                    col_name="ARO Accession",
                    samp_bin_name="sample1",
                ),
                call(
                    path=f"{tmp_dir}/sample1/output.gene_mapping_data.txt",
                    col_name="ARO Accession",
                    samp_bin_name="sample1",
                ),
                call(
                    path=f"{tmp_dir}/sample2/output.allele_mapping_data.txt",
                    col_name="ARO Accession",
                    samp_bin_name="sample2",
                ),
                call(
                    path=f"{tmp_dir}/sample2/output.gene_mapping_data.txt",
                    col_name="ARO Accession",
                    samp_bin_name="sample2",
                ),
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
        with patch("q2_amr.card.reads.run_command") as mock_run_command:
            run_rgi_bwt(
                "path_tmp",
                "sample1",
                "path_fwd",
                "path_rev",
                "bowtie2",
                8,
                True,
                True,
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
                    "--include_wildcard",
                    "--include_other_models",
                ],
                "path_tmp",
                verbose=True,
            )

    def test_exception_raised(self):
        expected_message = (
            "An error was encountered while running rgi, "
            "(return code 1), please inspect stdout and stderr to learn more."
        )

        with patch(
            "q2_amr.card.reads.run_command"
        ) as mock_run_command, self.assertRaises(Exception) as cm:
            mock_run_command.side_effect = subprocess.CalledProcessError(1, "cmd")
            run_rgi_bwt(
                cwd="path/cwd",
                samp="sample1",
                fwd="path/fwd",
                rev="path/rev",
                aligner="bwa",
                threads=1,
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
        with patch("q2_amr.card.reads.extract_sample_stats"), patch(
            "q2_amr.card.reads.plot_sample_stats",
            side_effect=mock_plot_sample_stats,
        ), tempfile.TemporaryDirectory() as tmp:
            visualize_annotation_stats(tmp, amr_reads_annotation)
            self.assertTrue(os.path.exists(os.path.join(tmp, "sample_stats_plot.html")))
            self.assertTrue(os.path.exists(os.path.join(tmp, "index.html")))
            self.assertTrue(os.path.exists(os.path.join(tmp, "q2templateassets")))
