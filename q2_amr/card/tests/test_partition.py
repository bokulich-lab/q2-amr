import os
import shutil

from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.partition import (
    collate_mags_annotations,
    collate_mags_kmer_analyses,
    collate_reads_allele_annotations,
    collate_reads_allele_kmer_analyses,
    collate_reads_gene_annotations,
    collate_reads_gene_kmer_analyses,
    partition_mags_annotations,
    partition_reads_allele_annotations,
    partition_reads_gene_annotations,
)
from q2_amr.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
)


class TestPartition(TestPluginBase):
    package = "q2_amr.card.tests"

    def test_partition_mags_annotations_duplicated(self):
        # Test partitioning mags annotations with duplicated mag ids
        self._test_partition_mags_annotations(
            duplicated=True, keys=[1, 2], bins=["1", "1"]
        )

    def test_partition_mags_annotations(self):
        # Test partitioning mags annotations with unique mag ids
        self._test_partition_mags_annotations(
            duplicated=False, keys=["bin2", "bin1"], bins=["2", "1"]
        )

    def test_partition_mags_warning_message(self):
        # Test warning message when partitioning mags annotations with num partitions
        # higher than the number of mags annotations
        annotations = self.setup_annotations(
            dir_name="collated/annotate_mags_output",
            format=CARDAnnotationDirectoryFormat,
        )
        with self.assertWarnsRegex(
            UserWarning, "You have requested a number of.*5.*2.*2"
        ):
            partition_mags_annotations(annotations=annotations, num_partitions=5)

    def _test_partition_mags_annotations(self, duplicated, keys, bins):
        # Set up for annotations
        annotations = self.setup_annotations(
            dir_name="collated/annotate_mags_output",
            format=CARDAnnotationDirectoryFormat,
        )
        # Change name of mag id so they are not duplicated
        if not duplicated:
            os.rename(
                os.path.join(annotations.path, "sample2", "bin1"),
                os.path.join(annotations.path, "sample2", "bin2"),
            )

        # Run partition_mags_annotations
        obs = partition_mags_annotations(annotations=annotations, num_partitions=2)

        # Assert if keys of collection are correct
        self.assertTrue(list(obs.keys()) == keys)

        # Assert if all files exist in the right location
        for key, samp, bin in zip(list(obs.keys()), ["2", "1"], bins):
            for file in ["amr_annotation.json", "amr_annotation.txt"]:
                path = os.path.join(obs[key].path, f"sample{samp}", f"bin{bin}", file)
                self.assertTrue(os.path.exists(path))

    def test_partition_reads_allele_annotations(self):
        self._test_partition_reads_annotations(
            dir="collated/annotate_reads_allele_output",
            files=["allele_mapping_data.txt", "overall_mapping_stats.txt"],
            format=CARDAlleleAnnotationDirectoryFormat,
            function=partition_reads_allele_annotations,
        )

    def test_partition_reads_gene_annotations(self):
        self._test_partition_reads_annotations(
            dir="collated/annotate_reads_gene_output",
            files=["gene_mapping_data.txt", "overall_mapping_stats.txt"],
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

    def test_collate_mags_annotations(self):
        # Test collate for mags annotations
        self._test_collate(
            data_dir="annotate_mags_output",
            files_to_assert=["bin1/amr_annotation.json", "bin1/amr_annotation.txt"],
            format=CARDAnnotationDirectoryFormat,
            function=collate_mags_annotations,
        )

    def test_collate_reads_allele_annotations(self):
        # Test collate for reads allele annotations
        self._test_collate(
            data_dir="annotate_reads_allele_output",
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
        self._test_collate(
            data_dir="annotate_reads_gene_output",
            files_to_assert=["gene_mapping_data.txt"],
            format=CARDGeneAnnotationDirectoryFormat,
            function=collate_reads_gene_annotations,
        )

    def test_collate_mags_kmer_analysis(self):
        # Test collate for MAGs k-mer analysis
        self._test_collate(
            data_dir="kmer_analysis_mags",
            files_to_assert=[
                "bin1/61mer_analysis.json",
                "bin1/61mer_analysis_rgi_summary.txt",
            ],
            format=CARDMAGsKmerAnalysisDirectoryFormat,
            function=collate_mags_kmer_analyses,
        )

    def test_collate_reads_allele_kmer_analysis(self):
        # Test collate for MAGs k-mer analysis
        self._test_collate(
            data_dir="kmer_analysis_reads_allele",
            files_to_assert=["61mer_analysis.json", "61mer_analysis.allele.txt"],
            format=CARDReadsAlleleKmerAnalysisDirectoryFormat,
            function=collate_reads_allele_kmer_analyses,
        )

    def test_collate_reads_gene_kmer_analysis(self):
        # Test collate for MAGs k-mer analysis
        self._test_collate(
            data_dir="kmer_analysis_reads_gene",
            files_to_assert=["61mer_analysis.json", "61mer_analysis.gene.txt"],
            format=CARDReadsGeneKmerAnalysisDirectoryFormat,
            function=collate_reads_gene_kmer_analyses,
        )

    def _test_collate(self, data_dir, files_to_assert, format, function):
        # Set up the list with annotations objects to collate
        artifact_1 = self.setup_annotations(
            dir_name=f"partitioned/{data_dir}_1", format=format
        )
        artifact_2 = self.setup_annotations(
            dir_name=f"partitioned/{data_dir}_2", format=format
        )

        artifacts = [artifact_1, artifact_2]

        # Run collate functions on the annotations
        collate = function(artifacts)

        # Assert if collated artifact has the correct format
        self.assertTrue(isinstance(collate, format))

        # Assert if all the files have been moved to the collated object
        for sample in ["sample1", "sample2"]:
            for file in files_to_assert:
                self.assertTrue(
                    os.path.exists(os.path.join(collate.path, sample, file))
                )

    def setup_annotations(self, dir_name, format):
        # Setup of the directory with dummy files and the needed directory format
        annotations = format()
        files = self.get_data_path(dir_name)
        shutil.copytree(files, annotations.path, dirs_exist_ok=True)
        return annotations
