import os
import shutil
import subprocess
import tarfile
from unittest.mock import MagicMock, call, patch

import requests
from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.database import fetch_card_db, preprocess
from q2_amr.types import CARDDatabaseDirectoryFormat, CARDKmerDatabaseDirectoryFormat


class TestAnnotateMagsCard(TestPluginBase):
    package = "q2_amr.tests"

    def mock_preprocess(self, dir, operation):
        if operation == "card":
            src_des = [
                ("DNA_fasta.fasta", "card_database_v3.2.7.fasta"),
                ("DNA_fasta_-.fasta", "card_database_v3.2.7_all.fasta"),
            ]
        else:
            src_des = [
                ("DNA_fasta.fasta", "wildcard_database_v0.fasta"),
                ("DNA_fasta_-.fasta", "wildcard_database_v0_all.fasta"),
            ]
        for file_name_src, file_name_des in src_des:
            shutil.copy(
                self.get_data_path(file_name_src), os.path.join(dir, file_name_des)
            )

    def test_fetch_card_db(self):
        f_card = open(self.get_data_path("card.tar.bz2"), "rb")
        f_wildcard = open(self.get_data_path("wildcard_data.tar.bz2"), "rb")
        mock_response_card = MagicMock(raw=f_card)
        mock_response_wildcard = MagicMock(raw=f_wildcard)

        with patch("requests.get") as mock_requests, patch(
            "q2_amr.card.database.preprocess", side_effect=self.mock_preprocess
        ):
            mock_requests.side_effect = [mock_response_card, mock_response_wildcard]

            obs = fetch_card_db()
        dir_file_list = [
            (str(obs[0]), "card.json"),
            (str(obs[0]), "card_database_v3.2.7.fasta"),
            (str(obs[0]), "card_database_v3.2.7_all.fasta"),
            (str(obs[0]), "wildcard_database_v0.fasta"),
            (str(obs[0]), "wildcard_database_v0_all.fasta"),
            (str(obs[0]), "index-for-model-sequences.txt"),
            (str(obs[1]), "61_kmer_db.json"),
            (str(obs[1]), "all_amr_61mers.txt"),
            (str(obs[0]), "nucleotide_fasta_protein_homolog_model_variants.fasta"),
            (
                str(obs[0]),
                "nucleotide_fasta_protein_overexpression_model_variants.fasta",
            ),
            (str(obs[0]), "nucleotide_fasta_protein_variant_model_variants.fasta"),
            (str(obs[0]), "nucleotide_fasta_rRNA_gene_variant_model_variants.fasta"),
        ]
        for dir, file in dir_file_list:
            self.assertTrue(os.path.exists(os.path.join(dir, file)))
        self.assertIsInstance(obs[0], CARDDatabaseDirectoryFormat)
        self.assertIsInstance(obs[1], CARDKmerDatabaseDirectoryFormat)

        expected_calls = [
            call("https://card.mcmaster.ca/latest/data", stream=True),
            call("https://card.mcmaster.ca/latest/variants", stream=True),
        ]
        mock_requests.assert_has_calls(expected_calls)

    def test_connection_error(self):
        with patch(
            "requests.get", side_effect=requests.ConnectionError
        ), self.assertRaisesRegex(
            requests.ConnectionError, "Network connectivity problems."
        ):
            fetch_card_db()

    def test_tarfile_read_error(self):
        with patch("tarfile.open", side_effect=tarfile.ReadError), patch(
            "requests.get"
        ), self.assertRaisesRegex(tarfile.ReadError, "Tarfile is invalid."):
            fetch_card_db()

    def test_subprocess_error(self):
        with patch(
            "q2_amr.card.database.run_command",
            side_effect=subprocess.CalledProcessError(1, "cmd"),
        ), self.assertRaisesRegex(
            Exception,
            "An error was encountered while running rgi, "
            r"\(return code 1\), please inspect stdout and stderr to learn more.",
        ):
            preprocess("path", "card")

    def test_preprocess_card(self):
        with patch("q2_amr.card.database.run_command") as mock_run_command:
            preprocess("path_tmp", "card")
            mock_run_command.assert_called_once_with(
                ["rgi", "card_annotation", "-i", "card/card.json"],
                "path_tmp",
                verbose=True,
            )

    def test_preprocess_wildcard(self):
        with patch("q2_amr.card.database.run_command") as mock_run_command:
            preprocess("path_tmp", "wildcard")
            mock_run_command.assert_called_once_with(
                [
                    "rgi",
                    "wildcard_annotation",
                    "-i",
                    "wildcard",
                    "--card_json",
                    "card/card.json",
                    "-v",
                    "0",
                ],
                "path_tmp",
                verbose=True,
            )
