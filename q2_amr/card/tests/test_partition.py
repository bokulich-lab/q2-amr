import os
import shutil

from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.partition import (
    partition_mags_annotations,
    partition_reads_allele_annotations,
    partition_reads_gene_annotations,
)
from q2_amr.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
)


class TestPartition(TestPluginBase):
    package = "q2_amr.card.tests"

    def test_partition_mags_annotations(self):
        # Set up for annotations
        annotations = self.setup_annotations(
            dir_name="collated/annotate_mags_output",
            format=CARDAnnotationDirectoryFormat,
        )

        # Run partition_mags_annotations
        obs = partition_mags_annotations(annotations=annotations, num_partitions=2)

        mag_ids = [
            "f5a16381-ea80-49f9-875e-620f333a9293",
            "e026af61-d911-4de3-a957-7e8bf837f30d",
        ]

        # Assert if keys of collection are correct
        self.assertTrue(list(obs.keys()) == mag_ids)

        # Assert if all files exist in the correct locations
        for mag_id, samp in zip(mag_ids, ["sample2", "sample1"]):
            for file in ["amr_annotation.json", "amr_annotation.txt"]:
                path = os.path.join(obs[mag_id].path, samp, mag_id, file)
                self.assertTrue(os.path.exists(path))

    def test_partition_mags_warning_message(self):
        # Test warning message when partitioning MAG annotations with num partitions
        # higher than the number of annotations
        annotations = self.setup_annotations(
            dir_name="collated/annotate_mags_output",
            format=CARDAnnotationDirectoryFormat,
        )
        with self.assertWarnsRegex(
            UserWarning, "You have requested a number of.*5.*2.*2"
        ):
            partition_mags_annotations(annotations=annotations, num_partitions=5)

    def setup_annotations(self, dir_name, format):
        # Setup of the directory with annotations files and the needed directory format
        annotations = format()
        files = self.get_data_path(dir_name)
        shutil.copytree(files, annotations.path, dirs_exist_ok=True)
        return annotations

    def test_partition_reads_allele_annotations(self):
        self._test_partition_reads_annotations(
            dir="collated/annotate_reads_allele_output",
            files=[
                "allele_mapping_data.txt",
                "overall_mapping_stats.txt",
                "sorted.length_100.bam",
            ],
            format=CARDAlleleAnnotationDirectoryFormat,
            function=partition_reads_allele_annotations,
        )

    def test_partition_reads_gene_annotations(self):
        self._test_partition_reads_annotations(
            dir="collated/annotate_reads_gene_output",
            files=["gene_mapping_data.txt"],
            format=CARDGeneAnnotationDirectoryFormat,
            function=partition_reads_gene_annotations,
        )

    def _test_partition_reads_annotations(self, dir, files, format, function):
        # Set up for annotations
        annotations = self.setup_annotations(dir, format)

        # Run function
        obs = function(annotations=annotations)

        # Assert if keys of collection are correct
        self.assertTrue(list(obs.keys()) == ["sample2", "sample1"])

        # Assert if all files exist in the right location
        for key, samp in zip(list(obs.keys()), ["2", "1"]):
            for file in files:
                path = os.path.join(obs[key].path, f"sample{samp}", file)
                self.assertTrue(os.path.exists(path))
