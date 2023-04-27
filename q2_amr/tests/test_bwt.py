import os
import tempfile

from q2_types.feature_data import DNAFASTAFormat
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase

from q2_amr.card import load_card_db, preprocess_card_db, load_card_db_fasta
from q2_amr.types import CARDDatabaseFormat


class TestBwt(TestPluginBase):
    package = 'q2_amr.tests'

    def test_load_card_db(self):
        card_db = Artifact.load(self.get_data_path('card_db_test.qza'))
        a = card_db.view(CARDDatabaseFormat)
        with tempfile.TemporaryDirectory() as tmp:
            load_card_db(tmp, a)
            assert os.path.exists(os.path.join(tmp, 'localDB', 'card.json'))
            assert os.path.exists(os.path.join(tmp, 'localDB', 'loaded_databases.json'))
            preprocess_card_db(tmp, a)
            assert os.path.exists(os.path.join(tmp, 'card_database_v3.2.5.fasta'))
            load_card_db_fasta(tmp, a)
            assert os.path.exists(os.path.join(tmp, 'localDB', 'card_reference.fasta'))

    def test_preprocess_card_db(self):
        card_db = Artifact.load(self.get_data_path('card_db.qza'))
        a = card_db.view(CARDDatabaseFormat)
        with tempfile.TemporaryDirectory() as tmp:
            preprocess_card_db(tmp, a)
            assert os.path.exists(os.path.join(tmp, 'card_database_v3.2.6.fasta'))
