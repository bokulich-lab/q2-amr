import os
import shutil
import subprocess
import tempfile
from unittest.mock import ANY, MagicMock, call, patch

from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from qiime2.plugin.testing import TestPluginBase

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
    CARDDatabaseDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
)


class TestAnnotateReadsCARD(TestPluginBase):
    package = "q2_amr.tests"

    @classmethod
    def setUpClass(cls):
        cls.sample_stats = {
            "sample1": {
                "total_reads": 5000,
                "mapped_reads": 59,
                "percentage": 1.18,
            },
            "sample2": {
                "total_reads": 7000,
                "mapped_reads": 212,
                "percentage": 3.03,
            },
        }

    def test_annotate_reads_card_single(self):
        self.annotate_reads_card_test_body("single")

    def test_annotate_reads_card_paired(self):
        self.annotate_reads_card_test_body("paired")

    def copy_needed_files(self, cwd, samp, **kwargs):
        # Create a sample directory
        samp_dir = os.path.join(cwd, samp)

        # Copy three dummy files to the directory
        for a, b in zip(["allele", "gene", "overall"], ["data", "data", "stats"]):
            shutil.copy(self.get_data_path(f"output.{a}_mapping_{b}.txt"), samp_dir)

    def annotate_reads_card_test_body(self, read_type):
        # Create single end or paired end reads object and CARD database object
        reads = (
            SingleLanePerSampleSingleEndFastqDirFmt()
            if read_type == "single"
            else SingleLanePerSamplePairedEndFastqDirFmt()
        )
        card_db = CARDDatabaseDirectoryFormat()

        # Copy manifest file to reads object
        manifest = self.get_data_path(f"MANIFEST_reads_{read_type}")
        shutil.copy(manifest, os.path.join(str(reads), "MANIFEST"))

        # Create MagicMock objects for run_rgi_bwt, run_rgi_load, read_in_txt and
        # create_count_table functions
        mock_run_rgi_bwt = MagicMock(side_effect=self.copy_needed_files)
        mock_run_rgi_load = MagicMock()
        mock_read_in_txt = MagicMock()
        mock_create_count_table = MagicMock()

        # Patch run_rgi_bwt, run_rgi_load, read_in_txt and create_count_table functions
        # and assign MagicMock objects
        with patch("q2_amr.card.reads.run_rgi_bwt", mock_run_rgi_bwt), patch(
            "q2_amr.card.reads.load_card_db", mock_run_rgi_load
        ), patch("q2_amr.card.reads.read_in_txt", mock_read_in_txt), patch(
            "q2_amr.card.reads.create_count_table", mock_create_count_table
        ):

            # Run annotate_reads_card function
            result = annotate_reads_card(reads, card_db)

            # Retrieve the path to cwd directory from mock_run_rgi_bwt arguments
            first_call_args = mock_run_rgi_bwt.call_args_list[0]
            tmp_dir = first_call_args.kwargs["cwd"]

            # Create four expected call objects for mock_run_rgi_bwt
            exp_calls_mock_bwt = [
                call(
                    cwd=tmp_dir,
                    aligner="kma",
                    threads=1,
                    include_wildcard=False,
                    include_other_models=False,
                    samp=f"sample{i}",
                    fwd=f"{reads}/sample{i}_00_L001_R1_001.fastq.gz",
                    rev=None
                    if read_type == "single"
                    else f"{reads}/sample{i}_00_L001_R2_001.fastq.gz",
                )
                for i in range(1, 3)
            ]

            # Expected call object for mock_run_rgi_load
            exp_calls_mock_load = [
                call(
                    tmp=tmp_dir,
                    card_db=ANY,
                    fasta=True,
                    include_other_models=False,
                    include_wildcard=False,
                ),
            ]

            # Create four expected call objects for mock_read_in_txt
            exp_calls_mock_read = [
                call(
                    path=f"{tmp_dir}/{samp}/output.{model}_mapping_data.txt",
                    samp_bin_name=samp,
                    data_type="reads",
                )
                for samp in ["sample1", "sample2"]
                for model in ["allele", "gene"]
            ]

            # Expected call objects for mock_create_count_table
            exp_calls_mock_count = [call([ANY, ANY]), call([ANY, ANY])]

            # Assert if all patched function were called with the expected calls
            mock_run_rgi_bwt.assert_has_calls(exp_calls_mock_bwt)
            mock_run_rgi_load.assert_has_calls(exp_calls_mock_load)
            mock_read_in_txt.assert_has_calls(exp_calls_mock_read)
            mock_create_count_table.assert_has_calls(exp_calls_mock_count)

            # Assert if all output files are the expected format
            self.assertIsInstance(result[0], CARDAlleleAnnotationDirectoryFormat)
            self.assertIsInstance(result[1], CARDGeneAnnotationDirectoryFormat)

            # Assert if the expected files are in every sample directory and in both
            # resulting CARD annotation objects
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
            run_rgi_bwt()
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

            self.assertEqual(sample_stats, self.sample_stats["sample1"])

    def test_plot_sample_stats(self):
        with tempfile.TemporaryDirectory() as tmp:
            plot_sample_stats(self.sample_stats, tmp)
            self.assertTrue(os.path.exists(os.path.join(tmp, "sample_stats_plot.html")))

    def mock_plot_sample_stats(self, sample_stats, output_dir):
        # Create a dummy HTML file and copy it to the output_dir
        with open(os.path.join(output_dir, "sample_stats_plot.html"), "w") as file:
            file.write("file")

    def test_visualize_annotation_stats(self):
        # Create a CARDGeneAnnotation object
        amr_reads_annotation = CARDGeneAnnotationDirectoryFormat()

        # Create two sample directories in the CARDGeneAnnotation object
        for num in range(1, 3):
            os.makedirs(os.path.join(str(amr_reads_annotation), f"sample{num}"))

        # Patch extract_sample_stats and plot_sample_stats with side effect
        # mock_plot_sample_stats
        with patch("q2_amr.card.reads.extract_sample_stats"), patch(
            "q2_amr.card.reads.plot_sample_stats",
            side_effect=self.mock_plot_sample_stats,
        ), tempfile.TemporaryDirectory() as tmp:

            # Run visualize_annotation_stats function
            visualize_annotation_stats(tmp, amr_reads_annotation)

            # Assert if all expected files are created
            for file in ["sample_stats_plot.html", "index.html", "q2templateassets"]:
                self.assertTrue(os.path.exists(os.path.join(tmp, file)))
