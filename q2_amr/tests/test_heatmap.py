import os
import tempfile
from unittest.mock import patch

from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.heatmap import change_names, run_rgi_heatmap


class TestHeatmap(TestPluginBase):
    package = "q2_amr.tests"

    # def test_heatmap(self):

    def test_run_rgi_heatmap(self):
        with patch("q2_amr.heatmap.run_command") as mock_run_command:
            run_rgi_heatmap(
                "path_tmp", "json_files_dir_path", "both", "drug_class", "fill", True
            )
            mock_run_command.assert_called_once_with(
                [
                    "rgi",
                    "heatmap",
                    "--input",
                    "json_files_dir_path",
                    "--output",
                    "path_tmp/results/heatmap",
                    "--display",
                    "fill",
                    "--clus",
                    "both",
                    "--cat",
                    "drug_class",
                    "--frequency",
                ],
                "path_tmp",
                verbose=True,
            )

    def test_change_names(self):
        extensions = [".eps", ".csv", ".png"]
        with tempfile.TemporaryDirectory() as tmp:
            results_dir = os.path.join(tmp, "results")
            os.makedirs(results_dir)
            for ext in extensions:
                with open(os.path.join(results_dir, f"heatmap22{ext}"), "w") as file:
                    file.write(f"{ext} file")
            change_names(results_dir)
            for ext in extensions:
                self.assertTrue(
                    os.path.exists(os.path.join(results_dir, f"heatmap{ext}"))
                )
