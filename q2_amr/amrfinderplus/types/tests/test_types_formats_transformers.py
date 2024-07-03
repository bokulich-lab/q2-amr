# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.core.exceptions import ValidationError
from qiime2.plugin.testing import TestPluginBase

from q2_amr.amrfinderplus.types._format import (
    AMRFinderPlusDatabaseDirectoryFormat,
    ARMFinderPlusAnnotationFormat,
)


class TestAMRFinderPlusTypesAndFormats(TestPluginBase):
    package = "q2_amr.amrfinderplus.types.tests"

    def test_amrfinderplus_database_directory_format_validate_positive(self):
        format = AMRFinderPlusDatabaseDirectoryFormat(
            self.get_data_path("database"), mode="r"
        )
        format.validate()

    def test_amrfinderplus_annotation_format_validate_positive(self):
        filepath = self.get_data_path("annotation/amr_annotation.tsv")
        format = ARMFinderPlusAnnotationFormat(filepath, mode="r")
        format.validate()

    def test_amrfinderplus_annotation_format_validate_positive_coordinates(self):
        filepath = self.get_data_path("annotation/amr_annotation_coordiantes.tsv")
        format = ARMFinderPlusAnnotationFormat(filepath, mode="r")
        format.validate()

    def test_amrfinderplus_annotation_format_validation_error(self):
        with self.assertRaises(ValidationError) as context:
            path = self.get_data_path("annotation/amr_annotation_wrong.tsv")
            format = ARMFinderPlusAnnotationFormat(path, mode="r")
            format.validate()

            header_coordinates = [
                "Protein identifier",
                "Contig id",
                "Start",
                "Stop",
                "Strand",
                "Gene symbol",
                "Sequence name",
                "Scope",
                "Element type",
                "Element subtype",
                "Class",
                "Subclass",
                "Method",
                "Target length",
                "Reference sequence length",
                "% Coverage of reference sequence",
                "% Identity to reference sequence",
                "Alignment length",
                "Accession of closest sequence",
                "Name of closest sequence",
                "HMM id",
                "HMM description",
            ]
            expected_message = (
                "Header line does not match ARMFinderPlusAnnotation format. Must "
                "consist of the following values: "
                + ", ".join(header_coordinates)
                + ".\nWhile Contig id, Start, Stop and Strand are optional."
                + "\n\nFound instead: "
                + "Incorrect Header 1, Incorrect Header 2, Incorrect Header 3"
            )

            self.assertEqual(str(context.exception), expected_message)
