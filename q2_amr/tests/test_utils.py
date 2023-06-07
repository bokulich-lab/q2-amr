from unittest.mock import patch

from qiime2.plugin.testing import TestPluginBase

from q2_amr.types import CARDDatabaseFormat
from q2_amr.utils import load_preprocess_card_db


class TestAnnotateReadsCARD(TestPluginBase):
    package = "q2_amr.tests"

    def test_load_card_db(self):
        card_db = CARDDatabaseFormat()
        with patch("q2_amr.utils.run_command") as mock_run_command:
            load_preprocess_card_db("path_tmp", card_db, "load")
            mock_run_command.assert_called_once_with(
                ["rgi", "load", "--card_json", str(card_db), "--local"],
                "path_tmp",
                verbose=True,
            )

    def test_preprocess_card_db(self):
        card_db = CARDDatabaseFormat()
        with patch("q2_amr.utils.run_command") as mock_run_command:
            load_preprocess_card_db("path_tmp", card_db, "preprocess")
            mock_run_command.assert_called_once_with(
                ["rgi", "card_annotation", "-i", str(card_db)], "path_tmp", verbose=True
            )

    def test_load_card_db_fasta(self):
        card_db = self.get_data_path("card_test.json")
        with patch("q2_amr.utils.run_command") as mock_run_command:
            load_preprocess_card_db("path_tmp", card_db, "load_fasta")
            mock_run_command.assert_called_once_with(
                [
                    "rgi",
                    "load",
                    "-i",
                    str(card_db),
                    "--card_annotation",
                    "card_database_v3.2.5.fasta",
                    "--local",
                ],
                "path_tmp",
                verbose=True,
            )
