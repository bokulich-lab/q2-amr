import os
import tarfile
from unittest.mock import MagicMock, patch

import requests
from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.database import fetch_card_db
from q2_amr.types import CARDDatabaseDirectoryFormat


class TestAnnotateMagsCard(TestPluginBase):
    package = "q2_amr.tests"

    def test_fetch_card_db(self):
        f = open(self.get_data_path("card.tar.bz2"), "rb")
        mock_response = MagicMock(raw=f)
        with patch("requests.get") as mock_requests:
            mock_requests.return_value = mock_response
            obs = fetch_card_db()
        self.assertTrue(os.path.exists(os.path.join(str(obs), "card.json")))
        self.assertIsInstance(obs, CARDDatabaseDirectoryFormat)
        mock_requests.assert_called_once_with(
            "https://card.mcmaster.ca/latest/data", stream=True
        )

    def test_fetch_card_data_connection_error(self):
        with patch(
            "requests.get", side_effect=requests.ConnectionError
        ), self.assertRaisesRegex(
            requests.ConnectionError, "Network connectivity problems."
        ):
            fetch_card_db()

    def test_fetch_card_data_tarfile_read_error(self):
        with patch(
            "tarfile.open", side_effect=tarfile.ReadError
        ), self.assertRaisesRegex(tarfile.ReadError, "Tarfile is invalid."):
            fetch_card_db()
