import tempfile

from qiime2.plugin.testing import TestPluginBase

from q2_amr.card import load_card_db


class TestBwt(TestPluginBase):
    package = 'q2_amr.types.tests'


    def test_load_card_db(self):
        card_database = self.get_data_path('card_db.qza')
        with tempfile.TemporaryDirectory() as tmp:
            load_card_db(tmp, card_database)
