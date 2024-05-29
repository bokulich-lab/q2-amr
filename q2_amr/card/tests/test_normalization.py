import os
import shutil
from unittest.mock import MagicMock, patch

import biom
from q2_types.feature_data import SequenceCharacteristicsDirectoryFormat
from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.normalization import normalize
from q2_amr.card.utils import InvalidParameterCombinationError


class TestNormalize(TestPluginBase):
    package = "q2_amr.card.tests"

    @classmethod
    def setUpClass(cls):
        # Mocking the biom.Table and gene_length class
        cls.table = MagicMock()
        cls.gene_length = MagicMock()

    def test_tpm_fpkm_uq_cuf_with_invalid_m_a_trim(self):
        # Test Error raised if gene-length is given with TMM method
        expected_message = (
            "Parameters m-trim and a-trim can only be used with methods TMM and " "CTF."
        )
        with self.assertRaises(InvalidParameterCombinationError) as cm:
            normalize(self.table, "tpm", m_trim=0.2, a_trim=0.05)
        self.assertEqual(str(cm.exception), expected_message)

    def test_tpm_fpkm_with_missing_gene_length(self):
        # Test Error raised if gene-length is missing with TPM method
        expected_message = "gene-length input is missing."
        with self.assertRaises(ValueError) as cm:
            normalize(self.table, "tpm")
        self.assertEqual(str(cm.exception), expected_message)

    def test_tmm_uq_cuf_ctf_with_gene_length(self):
        # Test Error raised if gene-length is given with TMM method
        expected_message = (
            "gene-length input can only be used with FPKM and TPM methods."
        )
        with self.assertRaises(ValueError) as cm:
            normalize(self.table, "tmm", gene_length=self.gene_length)
        self.assertEqual(str(cm.exception), expected_message)

    def test_tpm_fpkm_with_short_gene_length(self):
        # Test Error raised if gene-length is missing genes
        gene_length = SequenceCharacteristicsDirectoryFormat()
        shutil.copy(
            self.get_data_path("gene_length_short.txt"),
            os.path.join(gene_length.path, "gene_length.txt"),
        )
        table = biom.load_table(self.get_data_path("feature-table.biom"))
        expected_message = (
            "There are genes present in the FeatureTable that are not present "
            "in the gene-length input. Missing lengths for genes: "
            "{'ARO:3000027|ID:1757|Name:emrA|NCBI:AP009048.1'}"
        )
        with self.assertRaises(ValueError) as cm:
            normalize(table, "tpm", gene_length=gene_length)
        self.assertEqual(str(cm.exception), expected_message)

    @patch("q2_amr.card.normalization.TPM")
    def test_tpm_fpkm_with_valid_inputs(self, mock_tpm):
        # Test valid inputs for TPM method
        gene_length = SequenceCharacteristicsDirectoryFormat()
        shutil.copy(self.get_data_path("gene_length.txt"), gene_length.path)
        table = biom.load_table(self.get_data_path("feature-table.biom"))
        normalize(table=table, gene_length=gene_length, method="tpm")

    @patch("q2_amr.card.normalization.TMM")
    def test_tmm_uq_cuf_ctf_with_valid_inputs(self, mock_tmm):
        # Test valid inputs for TMM method
        table = biom.load_table(self.get_data_path("feature-table.biom"))
        normalize(table=table, method="tmm", a_trim=0.06, m_trim=0.4)
