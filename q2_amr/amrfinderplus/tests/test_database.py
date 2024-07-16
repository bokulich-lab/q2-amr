import os
import subprocess
from unittest.mock import patch

from qiime2.plugin.testing import TestPluginBase

from q2_amr.amrfinderplus.database import (
    _copy_all,
    fetch_amrfinderplus_db,
    run_amrfinder_u,
)


class TestFetchAMRFinderPlusDB(TestPluginBase):
    package = "q2_amr.amrfinderplus.tests"

    @patch("q2_amr.amrfinderplus.database.run_amrfinder_u")
    @patch("q2_amr.amrfinderplus.database._copy_all")
    def test_fetch_amrfinderplus_db(self, mock_run_amrfinder_u, mock__copy_all):
        fetch_amrfinderplus_db()

    @patch("q2_amr.amrfinderplus.database.run_command")
    def test_run_amrfinder_u(self, mock_run_command):
        run_amrfinder_u()
        mock_run_command.assert_called_once_with(
            ["amrfinder", "-u"],
            verbose=True,
        )

    @patch("q2_amr.amrfinderplus.database.run_command")
    def test_run_amrfinder_u_error(self, mock_run_command):
        expected_message = (
            "An error was encountered while running AMRFinderPlus, "
            "(return code 1), please inspect stdout and stderr to learn more."
        )
        mock_run_command.side_effect = subprocess.CalledProcessError(1, "cmd")
        with self.assertRaises(Exception) as cm:
            run_amrfinder_u()
        self.assertEqual(str(cm.exception), expected_message)

    def test__copy_all(self):
        tmp = self.temp_dir.name
        os.mkdir(os.path.join(tmp, "src"))
        os.mkdir(os.path.join(tmp, "des"))

        with open(os.path.join(tmp, "src", "a"), "w"), open(
            os.path.join(tmp, "src", "AMR_CDS.nto"), "w"
        ):
            pass

        _copy_all(os.path.join(tmp, "src"), os.path.join(tmp, "des"))
        self.assertTrue(os.path.exists(os.path.join(tmp, "des", "a")))
