import subprocess
from unittest.mock import patch

from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.utils import load_preprocess_card_db
from q2_amr.types import CARDDatabaseFormat


class TestAnnotateReadsCARD(TestPluginBase):
    package = "q2_amr.tests"

    def test_load_card_db(self):
        card_db = CARDDatabaseFormat()
        with patch("q2_amr.card.utils.run_command") as mock_run_command:
            load_preprocess_card_db("path_tmp", card_db, "load")
            mock_run_command.assert_called_once_with(
                ["rgi", "load", "--card_json", str(card_db), "--local"],
                "path_tmp",
                verbose=True,
            )

    def test_preprocess_card_db(self):
        card_db = CARDDatabaseFormat()
        with patch("q2_amr.card.utils.run_command") as mock_run_command:
            load_preprocess_card_db("path_tmp", card_db, "preprocess")
            mock_run_command.assert_called_once_with(
                ["rgi", "card_annotation", "-i", str(card_db)], "path_tmp", verbose=True
            )

    def test_load_card_db_fasta(self):
        card_db = self.get_data_path("card_test.json")
        with patch("q2_amr.card.utils.run_command") as mock_run_command:
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

    @patch("q2_amr.card.utils.run_command")
    def test_exception_raised(self, mock_run_command):
        mock_run_command.side_effect = subprocess.CalledProcessError(1, "cmd")
        tmp = "path/to/tmp"
        card_db = "path/to/card_db.json"
        expected_message = (
            "An error was encountered while running rgi, "
            "(return code 1), please inspect stdout and stderr to learn more."
        )
        operation = "load"
        with self.assertRaises(Exception) as cm:
            load_preprocess_card_db(tmp, card_db, operation)
        self.assertEqual(str(cm.exception), expected_message)
