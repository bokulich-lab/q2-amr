import os
import shutil
import tarfile
import tempfile

import requests

from q2_amr.types import CARDDatabaseDirectoryFormat

CARD_URL = "https://card.mcmaster.ca/latest/data"


def fetch_card_db() -> CARDDatabaseDirectoryFormat:
    try:
        response = requests.get(CARD_URL, stream=True)
    except requests.ConnectionError as e:
        raise requests.ConnectionError("Network connectivity problems.") from e
    with tempfile.TemporaryDirectory() as tmp_dir:
        try:
            with tarfile.open(fileobj=response.raw, mode="r|bz2") as tar:
                tar.extractall(path=tmp_dir)
        except tarfile.ReadError as a:
            raise tarfile.ReadError("Tarfile is invalid.") from a
        card_db = CARDDatabaseDirectoryFormat()
        shutil.move(
            os.path.join(tmp_dir, "card.json"), os.path.join(str(card_db), "card.json")
        )
        return card_db
