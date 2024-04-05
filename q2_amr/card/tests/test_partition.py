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

    def test_collate_mags_annotations(self):
        # Test collate for mags annotations
        self._test_collate(
            data_dir="annotate_mags_output",
            files_to_assert=["amr_annotation.json", "amr_annotation.txt"],
            samples=["sample1/bin1", "sample2/bin2"],
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
            samples=["sample1", "sample2"],
            format=CARDAlleleAnnotationDirectoryFormat,
            function=collate_reads_allele_annotations,
        )

    def test_collate_reads_gene_annotations(self):
        # Test collate for reads gene annotations
        self._test_collate(
            data_dir="annotate_reads_gene_output",
            files_to_assert=["gene_mapping_data.txt"],
            samples=["sample1", "sample2"],
            format=CARDGeneAnnotationDirectoryFormat,
            function=collate_reads_gene_annotations,
        )

    def test_collate_mags_kmer_analysis(self):
        # Test collate for MAGs k-mer analysis
        self._test_collate(
            data_dir="kmer_analysis_mags",
            files_to_assert=["61mer_analysis.json", "61mer_analysis_rgi_summary.txt"],
            samples=["sample1/bin1", "sample2/bin2"],
            format=CARDMAGsKmerAnalysisDirectoryFormat,
            function=collate_mags_kmer_analyses,
        )

    def test_collate_reads_allele_kmer_analysis(self):
        # Test collate for MAGs k-mer analysis
        self._test_collate(
            data_dir="kmer_analysis_reads_allele",
            files_to_assert=["61mer_analysis.json", "61mer_analysis.allele.txt"],
            samples=["sample1", "sample2"],
            format=CARDReadsAlleleKmerAnalysisDirectoryFormat,
            function=collate_reads_allele_kmer_analyses,
        )

    def test_collate_reads_gene_kmer_analysis(self):
        # Test collate for MAGs k-mer analysis
        self._test_collate(
            data_dir="kmer_analysis_reads_gene",
            files_to_assert=["61mer_analysis.json", "61mer_analysis.gene.txt"],
            samples=["sample1", "sample2"],
            format=CARDReadsGeneKmerAnalysisDirectoryFormat,
            function=collate_reads_gene_kmer_analyses,
        )

    def _test_collate(self, data_dir, files_to_assert, samples, format, function):
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
        for sample in samples:
            for file in files_to_assert:
                self.assertTrue(
                    os.path.exists(os.path.join(collate.path, sample, file))
                )

    def test_mags_file_exists_error(self):
        # Set up the list with duplicated artifacts
        artifact = self.setup_annotations(
            dir_name="partitioned/kmer_analysis_reads_allele_1",
            format=CARDReadsAlleleKmerAnalysisDirectoryFormat,
        )

        artifacts = [artifact, artifact]

        pattern = (
            r"The directory already exists: .*/sample1. Sample IDs must be "
            r"unique across all artifacts. Each artifact in the list must be "
            r"unique and cannot be repeated."
        )

        # Check if error is raised
        with self.assertRaisesRegex(FileExistsError, pattern):
            collate_reads_allele_kmer_analyses(artifacts)

    def test_reads_file_exists_error(self):
        # Set up the list with duplicated artifacts
        artifact = self.setup_annotations(
            dir_name="partitioned/annotate_mags_output_1",
            format=CARDAnnotationDirectoryFormat,
        )

        artifacts = [artifact, artifact]

        pattern = (
            r"The directory already exists: .*/bin1. MAG IDs must be "
            r"unique across all artifacts. Each artifact in the list must be "
            r"unique and cannot be repeated."
        )

        # Check if error is raised
        with self.assertRaisesRegex(FileExistsError, pattern):
            collate_reads_allele_kmer_analyses(artifacts)

    def setup_annotations(self, dir_name, format):
        # Setup of the directory with dummy files and the needed directory format
        annotations = format()
        files = self.get_data_path(dir_name)
        shutil.copytree(files, annotations.path, dirs_exist_ok=True)
        return annotations
