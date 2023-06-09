import os
import tarfile
from unittest.mock import MagicMock, patch

import requests
from qiime2.plugin.testing import TestPluginBase

from q2_amr.fetch_card_db import fetch_card_db
from q2_amr.types import CARDDatabaseDirectoryFormat


class TestAnnotateMagsCard(TestPluginBase):
    package = "q2_amr.tests"

    @patch("requests.get")
    def test_fetch_card_db(self, mock_requests):
        f = open(self.get_data_path("card.tar.bz2"), "rb")
        mock_response = MagicMock(raw=f)
        mock_requests.return_value = mock_response
        obs = fetch_card_db()
        self.assertTrue(os.path.exists(os.path.join(str(obs), "card.json")))
        self.assertIsInstance(obs, CARDDatabaseDirectoryFormat)
        mock_requests.assert_called_once_with(
            "https://card.mcmaster.ca/latest/data", stream=True
        )

    @patch("requests.get", side_effect=requests.ConnectionError)
    def test_fetch_card_data_connection_error(self, mock_requests):
        with self.assertRaisesRegex(
            requests.ConnectionError, "Network connectivity problems."
        ):
            fetch_card_db()

    @patch("tarfile.open", side_effect=tarfile.ReadError)
    def test_fetch_card_data_tarfile_read_error(self, mock_requests):
        with self.assertRaisesRegex(tarfile.ReadError, "Tarfile is invalid."):
            fetch_card_db()
