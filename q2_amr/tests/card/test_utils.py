import subprocess
from unittest.mock import patch

import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.utils import create_count_table, load_preprocess_card_db, read_in_txt
from q2_amr.types import CARDDatabaseFormat


class TestAnnotateReadsCARD(TestPluginBase):
    package = "q2_amr.tests"

    @classmethod
    def setUpClass(cls):
        cls.count_df_list = []
        for colname, ARG, sample in zip(
            ["ARO Term", "ARO Term", "Best_Hit_ARO"],
            ["mdtF", "mdtE", "mdtF"],
            ["sample1", "sample2", "sample1"],
        ):
            df = pd.DataFrame(
                {
                    colname: [ARG, "mgrA", "OprN", "mepA"],
                    sample: ["1", "1", "1", "1"],
                }
            )
            cls.count_df_list.append(df)

            cls.frequency_table = pd.DataFrame(
                {
                    "sample_id": ["sample1", "sample2"],
                    "mdtF": ["1", "0"],
                    "mgrA": ["1", "1"],
                    "OprN": ["1", "1"],
                    "mepA": ["1", "1"],
                    "mdtE": ["0", "1"],
                }
            )

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

    def test_exception_raised(self):
        tmp = "path/to/tmp"
        card_db = "path/to/card_db.json"
        expected_message = (
            "An error was encountered while running rgi, "
            "(return code 1), please inspect stdout and stderr to learn more."
        )
        operation = "load"
        with patch(
            "q2_amr.card.utils.run_command"
        ) as mock_run_command, self.assertRaises(Exception) as cm:
            mock_run_command.side_effect = subprocess.CalledProcessError(1, "cmd")
            load_preprocess_card_db(tmp, card_db, operation)
            self.assertEqual(str(cm.exception), expected_message)

    def test_read_in_txt_mags(self):
        # Test read_read_in_txt with output data from annotate_mags_card
        self.read_in_txt_test_body(
            "output.mags.txt", "sample1", self.count_df_list[2], "mags"
        )

    def test_read_in_txt_reads(self):
        # Test read_read_in_txt with output data from annotate_reads_card
        self.read_in_txt_test_body(
            "output.allele_mapping_data.txt", "sample1", self.count_df_list[0], "reads"
        )

    def read_in_txt_test_body(self, txt_file, samp_bin_name, mapping_data, data_type):
        # Create expected and observed count dataframes and compares them
        exp = mapping_data
        obs = read_in_txt(self.get_data_path(txt_file), samp_bin_name, data_type)
        pd.testing.assert_frame_equal(exp, obs)

    def test_create_count_table(self):
        # Create observed count table with create_count_table function
        df_list = [self.count_df_list[0], self.count_df_list[1]]
        obs = create_count_table(df_list)
        obs = obs.astype(str)

        # Create expected count table from
        exp = self.frequency_table
        exp.set_index("sample_id", inplace=True)

        # Compare expected and observed count table
        pd.testing.assert_frame_equal(exp, obs)

    def test_create_count_table_value_error(self):
        # Assert if ValueError is called when empy list is passed
        self.assertRaises(ValueError, create_count_table, [])
