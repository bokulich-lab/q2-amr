import os
import shutil
from unittest.mock import MagicMock, patch

import pandas as pd
from pandas._testing import assert_series_equal
from q2_types.feature_data import SequenceCharacteristicsDirectoryFormat
from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.normalization import _convert_lengths, _validate_parameters, normalize


class TestNormalize(TestPluginBase):
    package = "q2_amr.card.tests"

    @classmethod
    def setUpClass(cls):
        cls.gene_length = MagicMock()
        cls.lengths = pd.Series(
            {
                "ARO:3000026|ID:377|Name:mepA|NCBI:AY661734.1": 1356.0,
                "ARO:3000027|ID:1757|Name:emrA|NCBI:AP009048.1": 1173.0,
            },
            name="values",
        )
        cls.lengths.index.name = "index"

    def test_validate_parameters_uq_with_m_a_trim(self):
        # Test Error raised if gene-length is given with UQ method
        with self.assertRaisesRegex(
            ValueError,
            "Parameters m-trim and a-trim can only "
            "be used with methods TMM and CTF.",
        ):
            _validate_parameters("uq", 0.2, 0.05, None)

    def test_validate_parameters_tpm_missing_gene_length(self):
        # Test Error raised if gene-length is missing with TPM method
        with self.assertRaisesRegex(ValueError, "gene-length input is missing."):
            _validate_parameters("tpm", None, None, None)

    def test_validate_parameters_tmm_gene_length(self):
        # Test Error raised if gene-length is given with TMM method
        with self.assertRaisesRegex(
            ValueError, "gene-length input can only be used with FPKM and TPM methods."
        ):
            _validate_parameters("tmm", None, None, gene_length=self.gene_length)

    def test_validate_parameters_default_m_a_trim(self):
        # Test if m_trim and a_trim get set to default values if None
        m_trim, a_trim = _validate_parameters("tmm", None, None, None)
        self.assertEqual(m_trim, 0.3)
        self.assertEqual(a_trim, 0.05)

    def test_validate_parameters_m_a_trim(self):
        # Test if m_trim and a_trim are not modified if not None
        m_trim, a_trim = _validate_parameters("tmm", 0.1, 0.06, None)
        self.assertEqual(m_trim, 0.1)
        self.assertEqual(a_trim, 0.06)

    def test_convert_lengths_gene_length(self):
        # Test Error raised if gene-length is missing genes
        gene_length = SequenceCharacteristicsDirectoryFormat()
        shutil.copy(
            self.get_data_path("sequence_characteristics.tsv"), gene_length.path
        )
        table = pd.read_csv(
            self.get_data_path("feature-table.tsv"), sep="\t", index_col="ID"
        )

        obs = _convert_lengths(table, gene_length=gene_length)
        assert_series_equal(obs, self.lengths)

    def test_convert_lengths_short_gene_length(self):
        # Test Error raised if gene-length is missing genes
        gene_length = SequenceCharacteristicsDirectoryFormat()
        shutil.copy(
            self.get_data_path("sequence_characteristics_short.tsv"),
            os.path.join(gene_length.path, "sequence_characteristics.tsv"),
        )
        table = pd.read_csv(
            self.get_data_path("feature-table.tsv"), sep="\t", index_col="ID"
        )
        with self.assertRaisesRegex(
            ValueError,
            "There are genes present in the FeatureTable that are not present "
            "in the gene-length input. Missing lengths for genes: "
            "{'ARO:3000027|ID:1757|Name:emrA|NCBI:AP009048.1'}",
        ):
            _convert_lengths(table, gene_length=gene_length)

    @patch("q2_amr.card.normalization.TPM")
    def test_tpm_fpkm_with_valid_inputs(self, mock_tpm):
        # Test valid inputs for TPM method
        gene_length = SequenceCharacteristicsDirectoryFormat()
        shutil.copy(
            self.get_data_path("sequence_characteristics.tsv"), gene_length.path
        )
        table = pd.read_csv(
            self.get_data_path("feature-table.tsv"), sep="\t", index_col="ID"
        )
        normalize(table=table, gene_length=gene_length, method="tpm")

    @patch("q2_amr.card.normalization.TMM")
    def test_tmm_uq_cuf_ctf_with_valid_inputs(self, mock_tmm):
        # Test valid inputs for TMM method
        table = pd.read_csv(
            self.get_data_path("feature-table.tsv"), sep="\t", index_col="ID"
        )
        normalize(table=table, method="tmm", a_trim=0.06, m_trim=0.4)
