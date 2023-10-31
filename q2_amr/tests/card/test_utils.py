import os
import shutil
import subprocess
from unittest.mock import call, patch

import pandas as pd
from qiime2.plugin.testing import TestPluginBase
from test_mags import TestAnnotateMagsCard

from q2_amr.card.utils import create_count_table, load_card_db, read_in_txt
from q2_amr.types import CARDDatabaseDirectoryFormat, CARDKmerDatabaseDirectoryFormat


class TestAnnotateReadsCARD(TestPluginBase):
    package = "q2_amr.tests"

    mapping_data_sample1 = pd.DataFrame(
        {
            "ARO Accession": [3000796, 3000815, 3000805, 3000026],
            "sample1": [1, 1, 1, 1],
        }
    )

    mags_mapping_data_sample1 = pd.DataFrame(
        {
            "ARO": [3000796, 3000815, 3000805, 3000026],
            "sample1": [1, 1, 1, 1],
        }
    )

    mapping_data_sample2 = pd.DataFrame(
        {
            "ARO Accession": [3000797, 3000815, 3000805, 3000026],
            "sample2": [1, 1, 1, 2],
        }
    )

    def test_load_card_db_fasta(self):
        card_db = CARDDatabaseDirectoryFormat()
        kmer_db = CARDKmerDatabaseDirectoryFormat()
        card_json = self.get_data_path("card_test.json")
        shutil.copy(card_json, os.path.join(str(card_db), "card.json"))
        kmer_txt = self.get_data_path("kmer_txt_test.txt")
        shutil.copy(kmer_txt, os.path.join(str(kmer_db), "all_amr_61mers.txt"))
        kmer_json = self.get_data_path("kmer_json_test.json")
        shutil.copy(kmer_json, os.path.join(str(kmer_db), "61_kmer_db.json"))

        with patch("q2_amr.card.utils.run_command") as mock_run_command:
            load_card_db(
                tmp="path_tmp",
                card_db=card_db,
                kmer_db=kmer_db,
                kmer=True,
                fasta=True,
                include_other_models=False,
                include_wildcard=True,
            )
            load_card_db(
                tmp="path_tmp",
                card_db=card_db,
                kmer_db=kmer_db,
                kmer=True,
                fasta=True,
                include_other_models=True,
                include_wildcard=True,
            )
            expected_calls = [
                call(
                    [
                        "rgi",
                        "load",
                        "--card_json",
                        os.path.join(str(card_db), "card.json"),
                        "--local",
                        "--card_annotation",
                        os.path.join(str(card_db), "card_database_v3.2.5.fasta"),
                        "--wildcard_annotation",
                        os.path.join(str(card_db), "wildcard_database_v0.fasta"),
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
                ),
                call(
                    [
                        "rgi",
                        "load",
                        "--card_json",
                        os.path.join(str(card_db), "card.json"),
                        "--local",
                        "--card_annotation_all_models",
                        os.path.join(str(card_db), "card_database_v3.2.5_all.fasta"),
                        "--wildcard_annotation_all_models",
                        os.path.join(str(card_db), "wildcard_database_v0_all.fasta"),
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
                ),
            ]

            mock_run_command.assert_has_calls(expected_calls, any_order=False)

    def test_exception_raised(self):
        tmp = "path/to/tmp"
        card_db = "path/to/card_db.json"
        expected_message = (
            "An error was encountered while running rgi, "
            "(return code 1), please inspect stdout and stderr to learn more."
        )
        with patch(
            "q2_amr.card.utils.run_command"
        ) as mock_run_command, self.assertRaises(Exception) as cm:
            mock_run_command.side_effect = subprocess.CalledProcessError(1, "cmd")
            load_card_db(tmp, card_db)
            self.assertEqual(str(cm.exception), expected_message)

    def test_read_in_txt_mags(self):
        path = self.get_data_path("output.mags.txt")
        self.read_in_txt_test_body(
            path, "ARO", "sample1", self.mags_mapping_data_sample1
        )

    def test_read_in_txt_allele(self):
        path = self.get_data_path("output.allele_mapping_data.txt")
        self.read_in_txt_test_body(
            path, "ARO Accession", "sample1", self.mapping_data_sample1
        )

    def test_read_in_txt_gene(self):
        path = self.get_data_path("output.gene_mapping_data.txt")
        self.read_in_txt_test_body(
            path, "ARO Accession", "sample1", self.mapping_data_sample1
        )

    def read_in_txt_test_body(self, path, col_name, samp_bin_name, mapping_data):
        exp = mapping_data
        obs = read_in_txt(path, col_name, samp_bin_name)
        obs[col_name] = obs[col_name].astype(int)
        pd.testing.assert_frame_equal(exp, obs)

    def test_create_count_table(self):
        df_list = [self.mapping_data_sample1, self.mapping_data_sample2]
        obs = create_count_table(df_list)
        mag_test_class = TestAnnotateMagsCard()
        exp = mag_test_class.table
        exp.set_index("sample_id", inplace=True)
        exp = exp.astype(float)
        exp.columns = exp.columns.astype(float)
        pd.testing.assert_frame_equal(exp, obs)
        df_list_empty = []
        self.assertRaises(ValueError, create_count_table, df_list_empty)
