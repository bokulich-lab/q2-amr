import tempfile

from qiime2.plugin.testing import TestPluginBase

from q2_amr.card import load_card_db
import pkg_resources


class bwtTestPluginBase(TestPluginBase):
    package = 'q2_amr.types.tests'

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.TemporaryDirectory(prefix='q2-amr-test-temp-')

    def tearDown(self):
        self.temp_dir.cleanup()

    def get_data_path(self, filename):
        return pkg_resources.resource_filename(self.package, 'data/%s' % filename)


class TestCARDDatabaseTypesAndFormats(bwtTestPluginBase):

    def test_load_card_db(self):
        card_database = self.get_data_path('card_db.qza')
        with tempfile.TemporaryDirectory() as tmp:
            load_card_db(tmp, card_database)
