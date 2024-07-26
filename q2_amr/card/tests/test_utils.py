import os
import shutil
import subprocess
from unittest.mock import call, patch

from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.types import (
    CARDDatabaseDirectoryFormat,
    CARDKmerDatabaseDirectoryFormat,
)
from q2_amr.card.utils import copy_files, load_card_db


class TestAnnotateReadsCARD(TestPluginBase):
    package = "q2_amr.card.tests"

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
                kmer_size = load_card_db(
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
                    cmd=[
                        "rgi",
                        "load",
                        "--card_json",
                        os.path.join(str(card_db), "card.json"),
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
                    cwd=None,
                    verbose=True,
                )
                for flag, parameter in zip(flags, parameters)
            ]

            # Assert if the kmer_size is "61"
            assert kmer_size == "61"

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

    def test_copy_files(self):
        # Setup test files
        self.tmp = self.temp_dir.name

        file_path_1 = os.path.join(self.tmp, "DNA_fasta.fasta")
        file_path_2 = os.path.join(self.tmp, "DNA_fasta_-.fasta")

        shutil.copy(self.get_data_path("DNA_fasta.fasta"), self.tmp)
        shutil.copy(self.get_data_path("DNA_fasta_-.fasta"), self.tmp)

        # Call the function
        file_paths = [file_path_1, file_path_2]
        dst_path_components = [self.tmp, "dst_folder_1", "dst_folder_2"]

        copy_files(file_paths, *dst_path_components)

        # Assert if both files have been copied to the correct location
        dst_path_1 = os.path.join(
            self.tmp, "dst_folder_1", "dst_folder_2", "DNA_fasta.fasta"
        )
        dst_path_2 = os.path.join(
            self.tmp, "dst_folder_1", "dst_folder_2", "DNA_fasta_-.fasta"
        )

        self.assertTrue(os.path.exists(dst_path_1))
        self.assertTrue(os.path.exists(dst_path_2))
