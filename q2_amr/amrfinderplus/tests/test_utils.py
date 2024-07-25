from unittest.mock import patch

from qiime2.plugin.testing import TestPluginBase

from q2_amr.amrfinderplus.utils import run_amrfinderplus_n


class TestAnnotateMagsCard(TestPluginBase):
    package = "q2_amr.amrfinderplus.tests"

    @patch("q2_amr.amrfinderplus.utils.run_command")
    def test_run_amrfinderplus_n(self, mock_run_command):
        run_amrfinderplus_n(
            working_dir="path_dir",
            amrfinderplus_db="amrfinderplus_db",
            dna_sequences="dna_sequences",
            protein_sequences="protein_sequences",
            gff="gff",
            organism="Escherichia",
            plus=True,
            report_all_equal=True,
            ident_min=1,
            curated_ident=False,
            coverage_min=1,
            translation_table="11",
            threads=4,
        )
        mock_run_command.assert_called_once_with(
            [
                "amrfinder",
                "--database",
                "amrfinderplus_db",
                "-o",
                "path_dir/amr_annotations.tsv",
                "--print_node",
                "-n",
                "dna_sequences",
                "--nucleotide_output",
                "path_dir/amr_genes.fasta",
                "-p",
                "protein_sequences",
                "--protein_output",
                "path_dir/amr_proteins.fasta",
                "-g",
                "gff",
                "--threads",
                "4",
                "--organism",
                "Escherichia",
                "--mutation_all",
                "path_dir/amr_all_mutations.tsv",
                "--plus",
                "--report_all_equal",
                "--ident_min",
                "1",
                "--coverage_min",
                "1",
                "--translation_table",
                "11",
            ],
            "path_dir",
            verbose=True,
        )

    @patch("q2_amr.amrfinderplus.utils.run_command")
    def test_run_amrfinderplus_n_minimal(self, mock_run_command):
        run_amrfinderplus_n(
            working_dir="path_dir",
            amrfinderplus_db="amrfinderplus_db",
            dna_sequences=None,
            protein_sequences=None,
            gff=None,
            organism=None,
            plus=False,
            report_all_equal=False,
            ident_min=None,
            curated_ident=True,
            coverage_min=None,
            translation_table=None,
            threads=None,
        )
        mock_run_command.assert_called_once_with(
            [
                "amrfinder",
                "--database",
                "amrfinderplus_db",
                "-o",
                "path_dir/amr_annotations.tsv",
                "--print_node",
                "--ident_min",
                "-1",
            ],
            "path_dir",
            verbose=True,
        )
