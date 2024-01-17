import os
import shutil
import subprocess
from unittest.mock import call, patch

import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.utils import create_count_table, load_card_db, read_in_txt
from q2_amr.types import CARDDatabaseDirectoryFormat, CARDKmerDatabaseDirectoryFormat


class TestAnnotateReadsCARD(TestPluginBase):
    package = "q2_amr.card.tests"

    @classmethod
    def setUpClass(cls):
        cls.count_df_list = []
        for colname, values in zip(
            ["Reference Sequence", "ARO Term", "Best_Hit_ARO"],
            [
                [
                    "ARO:3000796|ID:121|Name:mdtF|NCBI:U00096.1",
                    "ARO:3000815|ID:154|Name:mgrA|NCBI:BA000018.3",
                    "ARO:3000805|ID:172|Name:OprN|NCBI:AE004091.2",
                    "ARO:3000026|ID:377|Name:mepA|NCBI:AY661734.1",
                ],
                ["mdtF", "mgrA", "OprN", "mepA"],
                ["mdtF", "mgrA", "OprN", "mepA"],
            ],
        ):
            df = pd.DataFrame(
                {
                    colname: values,
                    "sample1": ["1", "1", "1", "1"],
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
        cls.frequency_table.set_index("sample_id", inplace=True)

    def test_load_card_db_fasta(self):
        # Create CARD and Kmer database objects
        card_db = CARDDatabaseDirectoryFormat()
        kmer_db = CARDKmerDatabaseDirectoryFormat()

        # Tuples with source file name, destination file name and destination directory
        src_des_dir = [
            ("card_test.json", "card.json", card_db),
            ("kmer_txt_test.txt", "all_amr_61mers.txt", kmer_db),
            ("kmer_json_test.json", "61_kmer_db.json", kmer_db),
        ]

        # Copy files in src_des_dir to CARD and Kmer database objects
        for src, des, dir in src_des_dir:
            shutil.copy(self.get_data_path(src), os.path.join(str(dir), des))

        # Patch run_command
        with patch("q2_amr.card.utils.run_command") as mock_run_command:
            # Run load_card_db two times with include_other_models set to True and False
            for parameters in [False, True]:
                load_card_db(
                    tmp="path_tmp",
                    card_db=card_db,
                    kmer_db=kmer_db,
                    kmer=True,
                    fasta=True,
                    include_wildcard=True,
                    include_other_models=parameters,
                )

            # Create two expected call objects
            flags = ["", "_all_models"]
            parameters = ["", "_all"]

            expected_calls = [
                call(
                    [
                        "rgi",
                        "load",
                        "--card_json",
                        os.path.join(str(card_db), "card.json"),
                        "--local",
                        f"--card_annotation{flag}",
                        os.path.join(
                            str(card_db), f"card_database_v3.2.5{parameter}.fasta"
                        ),
                        f"--wildcard_annotation{flag}",
                        os.path.join(
                            str(card_db), f"wildcard_database_v0{parameter}.fasta"
                        ),
                        "--wildcard_index",
                        os.path.join(str(card_db), "index-for-model-sequences.txt"),
                        "--kmer_database",
                        os.path.join(str(kmer_db), "61_kmer_db.json"),
                        "--amr_kmers",
                        os.path.join(str(kmer_db), "all_amr_61mers.txt"),
                        "--kmer_size",
                        "61",
                    ],
                    "path_tmp",
                    verbose=True,
                )
                for flag, parameter in zip(flags, parameters)
            ]

            # Assert if function was called with expected calls
            mock_run_command.assert_has_calls(expected_calls, any_order=False)

    def test_exception_raised(self):
        # Simulate a subprocess.CalledProcessError during run_command
        expected_message = (
            "An error was encountered while running rgi, "
            "(return code 1), please inspect stdout and stderr to learn more."
        )
        with patch(
            "q2_amr.card.utils.run_command"
        ) as mock_run_command, self.assertRaises(Exception) as cm:
            mock_run_command.side_effect = subprocess.CalledProcessError(1, "cmd")
            load_card_db()
            self.assertEqual(str(cm.exception), expected_message)

    def test_read_in_txt_mags(self):
        # Test read_in_txt with output data from annotate_mags_card
        self.read_in_txt_test_body(
            filename="output.mags.txt",
            samp_bin_name="sample1",
            exp=self.count_df_list[2],
            data_type="mags",
        )

    def test_read_in_txt_reads_allele(self):
        # Test read_in_txt with allele mapping output data from annotate_reads_card
        self.read_in_txt_test_body(
            filename="output.allele_mapping_data.txt",
            samp_bin_name="sample1",
            exp=self.count_df_list[0],
            data_type="reads",
            map_type="allele",
        )

    def test_read_in_txt_reads_gene(self):
        # Test read_in_txt with gene mapping output data from annotate_reads_card
        self.read_in_txt_test_body(
            filename="output.gene_mapping_data.txt",
            samp_bin_name="sample1",
            exp=self.count_df_list[1],
            data_type="reads",
            map_type="gene",
        )

    def read_in_txt_test_body(
        self, filename, samp_bin_name, exp, data_type, map_type=None
    ):
        # Create expected and observed count dataframes and compare them
        obs = read_in_txt(
            self.get_data_path(filename), samp_bin_name, data_type, map_type
        )
        pd.testing.assert_frame_equal(exp, obs)

    def test_create_count_table(self):
        # Create list of dataframes to be used by create_count_table
        df_list = [self.count_df_list[1].copy(), self.count_df_list[1].copy()]
        df_list[1].iloc[0, 0] = "mdtE"
        df_list[1].rename(columns={"sample1": "sample2"}, inplace=True)

        # Create observed count table with create_count_table function
        obs = create_count_table(df_list)
        obs = obs.astype(str)

        # Define expected count table
        exp = self.frequency_table

        # Compare expected and observed count table
        pd.testing.assert_frame_equal(exp, obs)

    def test_create_count_table_value_error(self):
        # Assert if ValueError is called when empy list is passed
        self.assertRaises(ValueError, create_count_table, [])
