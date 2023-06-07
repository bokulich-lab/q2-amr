import os
import shutil
from unittest.mock import patch

from q2_types_genomics.per_sample_data import MultiMAGSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_amr.card import annotate_mags_card, run_rgi_main
from q2_amr.types import CARDAnnotationDirectoryFormat, CARDDatabaseFormat


class TestAnnotateMagsCard(TestPluginBase):
    package = "q2_amr.tests"

    def test_annotate_mags_card(self):
        output_txt = self.get_data_path("rgi_output.txt")
        output_json = self.get_data_path("rgi_output.json")
        manifest = self.get_data_path("MANIFEST_mags")
        mag = MultiMAGSequencesDirFmt()
        card_db = CARDDatabaseFormat()
        shutil.copy(manifest, os.path.join(str(mag), "MANIFEST"))

        def mock_run_rgi_main(
            tmp,
            input_sequence,
            alignment_tool,
            input_type,
            split_prodigal_jobs,
            include_loose,
            include_nudge,
            low_quality,
            num_threads,
        ):
            shutil.copy(output_txt, f"{tmp}/output.txt")
            shutil.copy(output_json, f"{tmp}/output.json")

        with patch("q2_amr.card.run_rgi_main", side_effect=mock_run_rgi_main), patch(
            "q2_amr.card.load_preprocess_card_db"
        ):
            result = annotate_mags_card(mag, card_db)
            self.assertIsInstance(result, CARDAnnotationDirectoryFormat)
            self.assertTrue(
                os.path.exists(
                    os.path.join(str(result), "sample1", "bin1", "amr_annotation.txt")
                )
            )
            self.assertTrue(
                os.path.exists(
                    os.path.join(str(result), "sample1", "bin1", "amr_annotation.json")
                )
            )

    def test_run_rgi_main(self):
        with patch("q2_amr.card.run_command") as mock_run_command:
            run_rgi_main(
                "path_tmp", "path_input", "DIAMOND", "contig", True, True, True, True, 8
            )
            mock_run_command.assert_called_once_with(
                [
                    "rgi",
                    "main",
                    "--input_sequence",
                    "path_input",
                    "--output_file",
                    "path_tmp/output",
                    "-n",
                    "8",
                    "--alignment_tool",
                    "DIAMOND",
                    "--input_type",
                    "contig",
                    "--local",
                    "--include_loose",
                    "--include_nudge",
                    "--low_quality",
                    "--split_prodigal_jobs",
                ],
                "path_tmp",
                verbose=True,
            )
