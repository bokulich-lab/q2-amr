import os
import shutil
from unittest.mock import MagicMock, patch

from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.kmer import (
    kmer_query_mags_card,
    kmer_query_reads_card,
    run_rgi_kmer_query,
)
from q2_amr.types import (
    CARDAnnotationDirectoryFormat,
    CARDDatabaseDirectoryFormat,
    CARDKmerDatabaseDirectoryFormat,
)
from q2_amr.types._format import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
)


class TestKmer(TestPluginBase):
    package = "q2_amr.card.tests"

    @classmethod
    def setUpClass(cls):
        cls.files_mags = (
            "output_61mer_analysis_rgi_summary.txt",
            "output_61mer_analysis.json",
        )

        cls.files_reads = (
            "output_61mer_analysis.allele.txt",
            "output_61mer_analysis.json",
            "output_61mer_analysis.gene.txt",
        )

    def copy_analysis_file(
        self, tmp, input_file, input_type, kmer_size, minimum, threads
    ):
        files = self.files_mags if input_type == "rgi" else self.files_reads
        for file in files:
            open(os.path.join(tmp, file), "w").close()

    def _run_kmer_query_test(self, annotation_format, output_format, query_function):
        mock_run_rgi_kmer_query = MagicMock(side_effect=self.copy_analysis_file)
        amr_annotations = annotation_format()
        card_db = CARDDatabaseDirectoryFormat()
        kmer_db = CARDKmerDatabaseDirectoryFormat()
        if query_function == kmer_query_reads_card:
            annotation_dir = os.path.join(amr_annotations.path, "sample1")
            os.makedirs(annotation_dir)
            des_path = os.path.join(annotation_dir, "sorted.length_100.bam")
            src_path = self.get_data_path("output.sorted.length_100.bam")
            shutil.copy(src_path, des_path)
            files = self.files_reads

        else:
            annotation_dir = os.path.join(amr_annotations.path, "sample1", "bin1")
            os.makedirs(annotation_dir)
            des_path = os.path.join(annotation_dir, "amr_annotation.json")
            src_path = self.get_data_path("rgi_output.json")
            shutil.copy(src_path, des_path)

            files = self.files_mags

        with patch(
            "q2_amr.card.kmer.run_rgi_kmer_query", side_effect=mock_run_rgi_kmer_query
        ), patch("q2_amr.card.kmer.load_card_db", return_value="61"):
            result = query_function(
                amr_annotations=amr_annotations, card_db=card_db, kmer_db=kmer_db
            )
            if query_function == kmer_query_reads_card:
                self.assertIsInstance(result[0], output_format[0])
                self.assertIsInstance(result[1], output_format[1])
                paths = [
                    os.path.join(str(result[0]), "sample1", files[0][7:]),
                    os.path.join(str(result[0]), "sample1", files[1][7:]),
                    os.path.join(str(result[1]), "sample1", files[2][7:]),
                ]
                for path in paths:
                    self.assertTrue(os.path.exists(path))
            else:
                self.assertIsInstance(result, output_format)

                for file in files:
                    path = os.path.join(str(result), "sample1", "bin1", file[7:])
                    self.assertTrue(os.path.exists(path))

    def test_kmer_query_mags_card(self):
        self._run_kmer_query_test(
            annotation_format=CARDAnnotationDirectoryFormat,
            output_format=CARDMAGsKmerAnalysisDirectoryFormat,
            query_function=kmer_query_mags_card,
        )

    def test_kmer_query_reads_card(self):
        self._run_kmer_query_test(
            annotation_format=CARDAlleleAnnotationDirectoryFormat,
            output_format=(
                CARDReadsAlleleKmerAnalysisDirectoryFormat,
                CARDReadsGeneKmerAnalysisDirectoryFormat,
            ),
            query_function=kmer_query_reads_card,
        )

    def test_run_rgi_kmer_query(self):
        with patch("q2_amr.card.kmer.run_command") as mock_run_command:
            run_rgi_kmer_query(
                tmp="path_tmp",
                input_file="path_input_file",
                input_type="rgi",
                kmer_size="61",
                minimum="10",
                threads="4",
            )
            mock_run_command.assert_called_once_with(
                [
                    "rgi",
                    "kmer_query",
                    "--input",
                    "path_input_file",
                    "--rgi",
                    "--kmer_size",
                    "61",
                    "--minimum",
                    "10",
                    "--threads",
                    "4",
                    "--output",
                    "output",
                    "--local",
                ],
                "path_tmp",
                verbose=True,
            )
