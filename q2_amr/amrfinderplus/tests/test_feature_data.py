import os
from pathlib import Path
from unittest.mock import MagicMock, patch

from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.genome_data import ProteinsDirectoryFormat
from qiime2.plugin.testing import TestPluginBase

from q2_amr.amrfinderplus.feature_data import (
    _get_file_paths,
    _move_or_create_files,
    _validate_inputs,
    annotate_feature_data_amrfinderplus,
)


class TestValidateInputs(TestPluginBase):
    package = "q2_amr.amrfinderplus.tests"

    def test_loci_mags(self):
        with self.assertRaisesRegex(
            ValueError,
            "Loci input can only be given in combination with proteins input.",
        ):
            _validate_inputs(mags="mags", loci="loci", proteins=None)

    def test_no_loci_protein_mags(self):
        with self.assertRaisesRegex(
            ValueError,
            "MAGs and proteins inputs together can only be given in combination with "
            "loci input.",
        ):
            _validate_inputs(mags="mags", loci=None, proteins="proteins")

    def test_no_protein_no_mags(self):
        with self.assertRaisesRegex(
            ValueError, "MAGs or proteins input has to be provided."
        ):
            _validate_inputs(mags=None, loci="loci_directory", proteins=None)


class TestMoveOrCreateFiles(TestPluginBase):
    package = "q2_amr.amrfinderplus.tests"

    def setUp(self):
        super().setUp()

        self.tmp = self.temp_dir.name
        self.src_dir = os.path.join(self.tmp, "src_dir")
        self.target_dir = os.path.join(self.tmp, "target_dir")
        os.mkdir(self.src_dir)
        os.mkdir(self.target_dir)

    def test_move_file(self):
        # Create a dummy file in the source directory
        with open(os.path.join(self.src_dir, "test_file.txt"), "w"):
            pass

        # Define the file operations
        file_operations = [("test_file.txt", self.target_dir)]

        # Run the function
        _move_or_create_files(
            src_dir=self.src_dir,
            mag_id="mag",
            file_operations=file_operations,
        )

        # Assert the file was moved
        self.assertTrue(
            os.path.exists(os.path.join(self.target_dir, "mag_test_file.txt"))
        )

    def test_file_missing_create_placeholder(self):
        # Define the file operations
        file_operations = [("test_file.txt", self.target_dir)]

        # Run the function
        _move_or_create_files(
            src_dir=self.src_dir,
            mag_id="mag",
            file_operations=file_operations,
        )

        # Assert the file was moved
        self.assertTrue(
            os.path.exists(os.path.join(self.target_dir, "mag_test_file.txt"))
        )

    def test_with_mags_and_proteins_file_missing(self):
        with self.assertRaisesRegex(
            ValueError, "Proteins file for ID 'mag_id' is missing in proteins input."
        ):
            _get_file_paths("path/mag_id.fasta", "path/mags", "path/proteins", None)


class TestGetFilePaths(TestPluginBase):
    package = "q2_amr.amrfinderplus.tests"

    def setUp(self):
        super().setUp()

        self.test_dir = self.temp_dir
        self.test_dir_path = Path(self.test_dir.name)
        self.file_path = self.test_dir_path / "test_file.fasta"
        self.file_path.touch()  # Create an empty test file

    def test_with_mags_and_proteins_file_exists(self):
        protein_file_path = self.test_dir_path / "test_file_proteins.fasta"
        protein_file_path.touch()  # Create an empty protein file

        mag_id, protein_path, gff_path = _get_file_paths(
            file=self.file_path,
            mags=self.test_dir_path,
            proteins=self.test_dir_path,
            loci=None,
        )
        self.assertEqual(mag_id, "test_file")
        self.assertEqual(protein_path, str(protein_file_path))
        self.assertIsNone(gff_path)

    def test_with_mags_and_proteins_file_missing(self):
        with self.assertRaisesRegex(
            ValueError,
            "Proteins file for ID 'test_file' is missing in proteins input.",
        ):
            _get_file_paths(
                file=self.file_path,
                mags=self.test_dir_path,
                proteins=self.test_dir_path,
                loci=None,
            )

    def test_with_proteins_only(self):
        protein_file_path = self.test_dir_path / "test_file_proteins.fasta"
        protein_file_path.touch()  # Create an empty protein file

        mag_id, protein_path, gff_path = _get_file_paths(
            file=protein_file_path, mags=None, proteins=self.test_dir_path, loci=None
        )
        self.assertEqual(mag_id, "test_file")
        self.assertEqual(protein_path, protein_file_path)
        self.assertIsNone(gff_path)

    def test_with_loci_file_exists(self):
        gff_file_path = self.test_dir_path / "test_file_loci.gff"
        gff_file_path.touch()  # Create an empty GFF file

        mag_id, protein_path, gff_path = _get_file_paths(
            file=self.file_path,
            mags=self.test_dir_path,
            proteins=None,
            loci=self.test_dir_path,
        )
        self.assertEqual(mag_id, "test_file")
        self.assertIsNone(protein_path)
        self.assertEqual(gff_path, str(gff_file_path))

    def test_with_loci_file_missing(self):
        with self.assertRaisesRegex(
            ValueError, "GFF file for ID 'test_file' is missing in loci input."
        ):
            _get_file_paths(
                file=self.file_path,
                mags=self.test_dir_path,
                proteins=None,
                loci=self.test_dir_path,
            )

    def test_with_mags_proteins_and_loci_all_files_exist(self):
        protein_file_path = self.test_dir_path / "test_file_proteins.fasta"
        gff_file_path = self.test_dir_path / "test_file_loci.gff"
        protein_file_path.touch()  # Create an empty protein file
        gff_file_path.touch()  # Create an empty GFF file

        mag_id, protein_path, gff_path = _get_file_paths(
            file=self.file_path,
            mags=self.test_dir_path,
            proteins=self.test_dir_path,
            loci=self.test_dir_path,
        )
        self.assertEqual(mag_id, "test_file")
        self.assertEqual(protein_path, str(protein_file_path))
        self.assertEqual(gff_path, str(gff_file_path))


class TestAnnotateFeatureDataAMRFinderPlus(TestPluginBase):
    package = "q2_amr.amrfinderplus.tests"

    @patch("q2_amr.amrfinderplus.feature_data._validate_inputs")
    @patch(
        "q2_amr.amrfinderplus.feature_data._get_file_paths",
        return_value=("mag_id", "protein_path", "gff_path"),
    )
    @patch("q2_amr.amrfinderplus.feature_data.run_amrfinderplus_n")
    @patch("q2_amr.amrfinderplus.feature_data._move_or_create_files")
    def test_annotate_feature_data_amrfinderplus_mags(
        self, mock_validate, mock_paths, mock_run, mock_move
    ):
        mags = MAGSequencesDirFmt()
        with open(os.path.join(str(mags), "mag.fasta"), "w"):
            pass
        annotate_feature_data_amrfinderplus(amrfinderplus_db=MagicMock(), mags=mags)

    @patch("q2_amr.amrfinderplus.feature_data._validate_inputs")
    @patch(
        "q2_amr.amrfinderplus.feature_data._get_file_paths",
        return_value=("mag_id", "protein_path", "gff_path"),
    )
    @patch("q2_amr.amrfinderplus.feature_data.run_amrfinderplus_n")
    @patch("q2_amr.amrfinderplus.feature_data._move_or_create_files")
    def test_annotate_feature_data_amrfinderplus_proteins(
        self, mock_validate, mock_paths, mock_run, mock_move
    ):
        proteins = ProteinsDirectoryFormat()
        with open(os.path.join(str(proteins), "proteins.fasta"), "w"):
            pass
        annotate_feature_data_amrfinderplus(
            amrfinderplus_db=MagicMock(), proteins=proteins
        )
