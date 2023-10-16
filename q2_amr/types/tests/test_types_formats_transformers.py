# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import tempfile

import pandas as pd
import pkg_resources
import qiime2
from Bio import SeqIO
from q2_types.feature_data import (
    DNAFASTAFormat,
    DNAIterator,
    ProteinFASTAFormat,
    ProteinIterator,
)
from q2_types_genomics.genome_data import GenesDirectoryFormat, ProteinsDirectoryFormat
from qiime2.plugin.testing import TestPluginBase
from skbio import DNA, Protein

from q2_amr.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
)
from q2_amr.types._format import (
    CARDAnnotationDirectoryFormat,
    CARDAnnotationTXTFormat,
    CARDDatabaseFormat,
)
from q2_amr.types._transformer import (
    _read_from_card_file,
    card_annotation_df_to_fasta,
    extract_sequence,
)


class AMRTypesTestPluginBase(TestPluginBase):
    package = "q2_amr.types.tests"

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.TemporaryDirectory(prefix="q2-amr-test-temp-")

    def tearDown(self):
        self.temp_dir.cleanup()

    def get_data_path(self, filename):
        return pkg_resources.resource_filename(self.package, "data/%s" % filename)


class TestCARDDatabaseTypesAndFormats(AMRTypesTestPluginBase):
    def test_card_database_format_validate_positive(self):
        filepath = self.get_data_path("card_test.json")
        format = CARDDatabaseFormat(filepath, mode="r")
        format.validate()

    def test_dataframe_to_card_format_transformer(self):
        filepath = self.get_data_path("card_test.json")
        transformer = self.get_transformer(pd.DataFrame, CARDDatabaseFormat)
        card_df = pd.read_json(filepath).transpose()
        obs = transformer(card_df)
        self.assertIsInstance(obs, CARDDatabaseFormat)

    def test_card_format_to_dataframe_transformer(self):
        filepath = self.get_data_path("card_test.json")
        transformer = self.get_transformer(CARDDatabaseFormat, pd.DataFrame)
        card_db = CARDDatabaseFormat(filepath, mode="r")
        obs = transformer(card_db)
        self.assertIsInstance(obs, pd.DataFrame)

    def test_protein_card_iterator_transformer(self):
        filepath = self.get_data_path("card_test.json")
        transformer = self.get_transformer(CARDDatabaseFormat, ProteinIterator)
        card_db = CARDDatabaseFormat(filepath, mode="r")
        obs = transformer(card_db)
        self.assertIsInstance(obs, ProteinIterator)

    def test_dna_card_iterator_transformer(self):
        filepath = self.get_data_path("card_test.json")
        transformer = self.get_transformer(CARDDatabaseFormat, DNAIterator)
        card_db = CARDDatabaseFormat(filepath, mode="r")
        obs = transformer(card_db)
        test_fasta = DNAFASTAFormat(self.get_data_path("card_test_dna.fasta"), mode="r")
        exp = test_fasta.view(DNAIterator)
        self.assertIsInstance(obs, DNAIterator)
        for e, o in zip(exp, obs):
            self.assertEqual(e, o)

    def test_protein_card_fasta_transformer(self):
        filepath = self.get_data_path("card_test.json")
        filepath2 = self.get_data_path("card_test_protein.fasta")
        transformer = self.get_transformer(CARDDatabaseFormat, ProteinFASTAFormat)
        card_db = CARDDatabaseFormat(filepath, mode="r")
        obs = transformer(card_db)
        obs_dict = {rec.id: rec.seq for rec in SeqIO.parse(str(obs), "fasta")}
        exp_dict = {rec.id: rec.seq for rec in SeqIO.parse(filepath2, "fasta")}
        self.assertDictEqual(obs_dict, exp_dict)
        self.assertIsInstance(obs, ProteinFASTAFormat)

    def test_dna_card_fasta_transformer(self):
        filepath = self.get_data_path("card_test.json")
        filepath2 = self.get_data_path("card_test_dna.fasta")
        transformer = self.get_transformer(CARDDatabaseFormat, DNAFASTAFormat)
        card_db = CARDDatabaseFormat(filepath, mode="r")
        obs = transformer(card_db)
        obs_dict = {rec.id: rec.seq for rec in SeqIO.parse(str(obs), "fasta")}
        exp_dict = {rec.id: rec.seq for rec in SeqIO.parse(filepath2, "fasta")}
        self.assertDictEqual(obs_dict, exp_dict)
        self.assertIsInstance(obs, DNAFASTAFormat)

    def test_extract_sequence_dna(self):
        with open(self.get_data_path("card_test.json"), "rb") as f:
            db = json.load(f)
        obs = extract_sequence("dna", "2", "1188", db)
        exp = DNA(
            "ATGAAAGCATATTTCATCGCCATACTTACCTTATTCACTTGTATAGCTACCGTCGTCCGGGCGCAGCAAATGTC"
            "TGAACTTGAAAACCGGATTGACAGTCTGCTCAATGGCAAGAAAGCCACCGTTGGTATAGCCGTATGGACAGACA"
            "AAGGAGACATGCTCCGGTATAACGACCATGTACACTTCCCCTTGCTCAGTGTATTCAAATTCCATGTGGCACTG"
            "GCCGTACTGGACAAGATGGATAAGCAAAGCATCAGTCTGGACAGCATTGTTTCCATAAAGGCATCCCAAATGCC"
            "GCCCAATACCTACAGCCCCCTGCGGAAGAAGTTTCCCGACCAGGATTTCACGATTACGCTTAGGGAACTGATGC"
            "AATACAGCATTTCCCAAAGCGACAACAATGCCTGCGACATCTTGATAGAATATGCAGGAGGCATCAAACATATC"
            "AACGACTATATCCACCGGTTGAGTATCGACTCCTTCAACCTCTCGGAAACAGAAGACGGCATGCACTCCAGCTT"
            "CGAGGCTGTATACCGCAACTGGAGTACTCCTTCCGCTATGGTCCGACTACTGAGAACGGCTGATGAAAAAGAGT"
            "TGTTCTCCAACAAGGAGCTGAAAGACTTCTTGTGGCAGACCATGATAGATACTGAAACCGGTGCCAACAAACTG"
            "AAAGGTATGTTGCCAGCCAAAACCGTGGTAGGACACAAGACCGGCTCTTCCGACCGCAATGCCGACGGTATGAA"
            "AACTGCAGATAATGATGCCGGCCTCGTTATCCTTCCCGACGGCCGGAAATACTACATTGCCGCCTTCGTCATGG"
            "ACTCATACGAGACGGATGAGGACAATGCGAACATCATCGCCCGCATATCACGCATGGTATATGATGCGATGAGA"
            "TGA"
        )
        exp.metadata["id"] = "gb|GQ343019.1|+|132-1023|ARO:3002999|CblA-1"
        exp.metadata["description"] = "[mixed culture bacterium AX_gF3SD01_15]"
        self.assertEqual(exp, obs)
        self.assertIsInstance(obs, DNA)

    def test_extract_sequence_protein(self):
        with open(self.get_data_path("card_test.json"), "rb") as f:
            db = json.load(f)
        obs = extract_sequence("protein", "2", "1188", db)
        exp = Protein(
            "MKAYFIAILTLFTCIATVVRAQQMSELENRIDSLLNGKKATVGIAVWTDKGDMLRYNDHVHFPLLSVFKFHVAL"
            "AVLDKMDKQSISLDSIVSIKASQMPPNTYSPLRKKFPDQDFTITLRELMQYSISQSDNNACDILIEYAGGIKHI"
            "NDYIHRLSIDSFNLSETEDGMHSSFEAVYRNWSTPSAMVRLLRTADEKELFSNKELKDFLWQTMIDTETGANKL"
            "KGMLPAKTVVGHKTGSSDRNADGMKTADNDAGLVILPDGRKYYIAAFVMDSYETDEDNANIIARISRMVYDAMR"
        )
        exp.metadata["id"] = "gb|ACT97415.1|ARO:3002999|CblA-1"
        exp.metadata["description"] = "[mixed culture bacterium AX_gF3SD01_15]"
        self.assertEqual(exp, obs)
        self.assertIsInstance(obs, Protein)

    def test_read_from_card_file(self):
        path = self.get_data_path("card_test.json")
        genomes = _read_from_card_file(path, "protein")
        self.assertIsInstance(genomes, ProteinFASTAFormat)

    def test_read_from_card_generator(self):
        path = self.get_data_path("card_test.json")
        genomes = _read_from_card_file(path, "protein")
        generator = ProteinIterator(genomes)
        self.assertIsInstance(generator, ProteinIterator)


class TestCARDMagsAnnotationTypesAndFormats(AMRTypesTestPluginBase):
    def test_df_to_card_annotation_format_transformer(self):
        filepath = self.get_data_path("rgi_output.txt")
        transformer = self.get_transformer(pd.DataFrame, CARDAnnotationTXTFormat)
        df = pd.read_csv(filepath, sep="\t")
        obs = transformer(df)
        self.assertIsInstance(obs, CARDAnnotationTXTFormat)

    def test_card_annotation_format_to_df_transformer(self):
        filepath = self.get_data_path("rgi_output.txt")
        transformer = self.get_transformer(CARDAnnotationTXTFormat, pd.DataFrame)
        card_anno = CARDAnnotationTXTFormat(filepath, mode="r")
        obs = transformer(card_anno)
        self.assertIsInstance(obs, pd.DataFrame)

    def test_card_annotation_df_to_fasta(self):
        filepath = self.get_data_path("rgi_output.txt")
        filepath2 = self.get_data_path("rgi_output_protein.fna")
        filepath3 = self.get_data_path("rgi_output_dna.fna")
        protein = card_annotation_df_to_fasta(filepath, "Protein")
        dna = card_annotation_df_to_fasta(filepath, "DNA")
        with open(str(protein), "r") as protein_fh_obs:
            protein_contents_obs = protein_fh_obs.read()
        with open(str(dna), "r") as dna_fh_obs:
            dna_contents_obs = dna_fh_obs.read()
        with open(filepath2, "r") as protein_fh_exp:
            protein_contents_exp = protein_fh_exp.read()
        with open(filepath3, "r") as dna_fh_exp:
            dna_contents_exp = dna_fh_exp.read()
        self.assertEqual(protein_contents_obs, protein_contents_exp)
        self.assertEqual(dna_contents_obs, dna_contents_exp)

    def test_CARDAnnotationDirectoryFormat_to_GenesDirectoryFormat_transformer(self):
        filepath = self.get_data_path("annotate_mags_output")
        transformer = self.get_transformer(
            CARDAnnotationDirectoryFormat, GenesDirectoryFormat
        )
        obs = transformer(filepath)
        self.assertIsInstance(obs, GenesDirectoryFormat)
        self.assertTrue(
            os.path.exists(os.path.join(str(obs), "sample1", "bin1_genes.fasta"))
        )
        self.assertTrue(
            os.path.exists(os.path.join(str(obs), "sample2", "bin1_genes.fasta"))
        )

    def test_CARDAnnotationDirectoryFormat_to_ProteinsDirectoryFormat_transformer(self):
        filepath = self.get_data_path("annotate_mags_output")
        transformer = self.get_transformer(
            CARDAnnotationDirectoryFormat, ProteinsDirectoryFormat
        )
        obs = transformer(filepath)
        self.assertIsInstance(obs, ProteinsDirectoryFormat)
        self.assertTrue(
            os.path.exists(os.path.join(str(obs), "sample1", "bin1_proteins.fasta"))
        )
        self.assertTrue(
            os.path.exists(os.path.join(str(obs), "sample2", "bin1_proteins.fasta"))
        )

    def test_CARDAnnotationDirectoryFormat_to_qiime2_Metadata_transformer(self):
        transformer = self.get_transformer(
            CARDAnnotationDirectoryFormat, qiime2.Metadata
        )
        annotation = CARDAnnotationDirectoryFormat(
            self.get_data_path("annotate_mags_output"), "r"
        )
        metadata_obt = transformer(annotation)
        self.assertIsInstance(metadata_obt, qiime2.Metadata)


class TestCARDReadsAnnotationTypesAndFormats(AMRTypesTestPluginBase):
    def test_CARDGeneAnnotationDirectoryFormat_to_qiime2_Metadata_transformer(self):
        transformer = self.get_transformer(
            CARDGeneAnnotationDirectoryFormat, qiime2.Metadata
        )
        annotation = CARDGeneAnnotationDirectoryFormat(
            self.get_data_path("annotate_reads_output"), "r"
        )
        metadata_obt = transformer(annotation)
        self.assertIsInstance(metadata_obt, qiime2.Metadata)

    def test_CARDAlleleAnnotationDirectoryFormat_to_qiime2_Metadata_transformer(self):
        transformer = self.get_transformer(
            CARDGeneAnnotationDirectoryFormat, qiime2.Metadata
        )
        annotation = CARDAlleleAnnotationDirectoryFormat(
            self.get_data_path("annotate_reads_output"), "r"
        )
        metadata_obt = transformer(annotation)
        self.assertIsInstance(metadata_obt, qiime2.Metadata)
