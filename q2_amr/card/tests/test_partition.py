import os
import shutil

from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.partition import (
    collate_mags_annotations,
    collate_reads_allele_annotations,
    collate_reads_gene_annotations,
)
from q2_amr.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
)


class TestPartition(TestPluginBase):
    package = "q2_amr.card.tests"

    def test_collate_mags_annotations(self):
        # Test collate for mags annotations
        self._test_collate_mags_annotations(
            annotations_dir="annotate_mags_output",
            files_to_assert=["bin1/amr_annotation.json", "bin1/amr_annotation.txt"],
            format=CARDAnnotationDirectoryFormat,
            function=collate_mags_annotations,
        )

    def test_collate_reads_allele_annotations(self):
        # Test collate for reads allele annotations
        self._test_collate_mags_annotations(
            annotations_dir="annotate_reads_allele_output",
            files_to_assert=[
                "allele_mapping_data.txt",
                "overall_mapping_stats.txt",
                "sorted.length_100.bam",
            ],
            format=CARDAlleleAnnotationDirectoryFormat,
            function=collate_reads_allele_annotations,
        )

    def test_collate_reads_gene_annotations(self):
        # Test collate for reads gene annotations
        self._test_collate_mags_annotations(
            annotations_dir="annotate_reads_gene_output",
            files_to_assert=["gene_mapping_data.txt"],
            format=CARDGeneAnnotationDirectoryFormat,
            function=collate_reads_gene_annotations,
        )

    def _test_collate_mags_annotations(
        self, annotations_dir, files_to_assert, format, function
    ):
        # Set up the list with annotations objects to collate
        annotation_1 = self.setup_annotations(
            dir_name=f"partitioned/{annotations_dir}_1", format=format
        )
        annotation_2 = self.setup_annotations(
            dir_name=f"partitioned/{annotations_dir}_2", format=format
        )

        annotations = [annotation_1, annotation_2]

        # Run collate functions on the annotations
        collate = function(annotations)

        # Assert if all the files have been moved to the collated object
        for sample in ["sample1", "sample2"]:
            for file in files_to_assert:
                self.assertTrue(
                    os.path.exists(os.path.join(collate.path, sample, file))
                )

    def setup_annotations(self, dir_name, format):
        # Setup of the directory with annotations files and the needed directory format
        annotations = format()
        files = self.get_data_path(dir_name)
        shutil.copytree(files, annotations.path, dirs_exist_ok=True)
        return annotations
