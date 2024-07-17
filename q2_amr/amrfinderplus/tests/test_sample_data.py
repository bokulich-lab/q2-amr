import os
from unittest.mock import MagicMock, patch

from q2_types.per_sample_sequences import ContigSequencesDirFmt, MultiMAGSequencesDirFmt
from qiime2.plugin.testing import TestPluginBase

from q2_amr.amrfinderplus.sample_data import annotate_sample_data_amrfinderplus
from q2_amr.amrfinderplus.types import AMRFinderPlusDatabaseDirFmt


def mock_run_amrfinderplus_n(
    working_dir,
    amrfinderplus_db,
    dna_sequence,
    protein_sequence,
    gff,
    organism,
    plus,
    report_all_equal,
    ident_min,
    curated_ident,
    coverage_min,
    translation_table,
    threads,
):
    with open(os.path.join(working_dir, "amr_annotations.tsv"), "w"):
        pass
    if organism:
        with open(os.path.join(working_dir, "amr_mutations.tsv"), "w"):
            pass
    if dna_sequence:
        with open(os.path.join(working_dir, "amr_genes.fasta"), "w"):
            pass
    if protein_sequence:
        with open(os.path.join(working_dir, "amr_proteins.fasta"), "w"):
            pass


class TestAnnotateSampleDataAMRFinderPlus(TestPluginBase):
    package = "q2_amr.amrfinderplus.tests"

    files_contigs = [
        "amr_annotations.tsv",
        "amr_mutations.tsv",
        "sample1_amr_genes.fasta",
    ]

    files_mags = [
        "mag1_amr_annotations.tsv",
        "mag1_amr_mutations.tsv",
        "mag1_amr_genes.fasta",
    ]

    def test_annotate_sample_data_amrfinderplus_mags(self):
        sequences = MultiMAGSequencesDirFmt()
        with open(os.path.join(str(sequences), "MANIFEST"), "w") as file:
            file.write("sample-id,mag-id,filename\nsample1,mag1,sample1/mag1.fasta\n")
        self._helper(sequences=sequences, organism=None, files=self.files_mags)

    def test_annotate_sample_data_amrfinderplus_mags_organism(self):
        sequences = MultiMAGSequencesDirFmt()
        with open(os.path.join(str(sequences), "MANIFEST"), "w") as file:
            file.write("sample-id,mag-id,filename\nsample1,mag1,sample1/mag1.fasta\n")
        self._helper(sequences, "Escherichia", files=self.files_mags)

    def test_annotate_sample_data_amrfinderplus_contigs(self):
        sequences = ContigSequencesDirFmt()
        with open(os.path.join(str(sequences), "sample1_contigs.fasta"), "w"):
            pass
        self._helper(sequences=sequences, organism=None, files=self.files_contigs)

    def test_annotate_sample_data_amrfinderplus_contigs_organism(self):
        sequences = ContigSequencesDirFmt()
        with open(os.path.join(str(sequences), "sample1_contigs.fasta"), "w"):
            pass
        self._helper(
            sequences=sequences, organism="Escherichia", files=self.files_contigs
        )

    def _helper(self, sequences, organism, files):
        amrfinderplus_db = AMRFinderPlusDatabaseDirFmt()
        mock_create_count_table = MagicMock()
        mock_read_in_txt = MagicMock()
        with patch(
            "q2_amr.amrfinderplus.sample_data.run_amrfinderplus_n",
            side_effect=mock_run_amrfinderplus_n,
        ), patch(
            "q2_amr.amrfinderplus.sample_data.read_in_txt", mock_read_in_txt
        ), patch(
            "q2_amr.amrfinderplus.sample_data.create_count_table",
            mock_create_count_table,
        ):
            result = annotate_sample_data_amrfinderplus(
                sequences=sequences,
                amrfinderplus_db=amrfinderplus_db,
                organism=organism,
            )
            self.assertTrue(
                os.path.exists(os.path.join(str(result[0]), "sample1", files[0]))
            )
            self.assertTrue(
                os.path.exists(os.path.join(str(result[1]), "sample1", files[1]))
            )
            self.assertTrue(os.path.exists(os.path.join(str(result[2]), files[2])))
