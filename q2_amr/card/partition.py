import itertools
import os
import warnings
from typing import Union

from qiime2.util import duplicate

from q2_amr.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
)


def collate_mags_annotations(
    annotations: CARDAnnotationDirectoryFormat,
) -> CARDAnnotationDirectoryFormat:
    return _collate(annotations)


def collate_reads_allele_annotations(
    annotations: CARDAlleleAnnotationDirectoryFormat,
) -> CARDAlleleAnnotationDirectoryFormat:
    return _collate(annotations)


def collate_reads_gene_annotations(
    annotations: CARDGeneAnnotationDirectoryFormat,
) -> CARDGeneAnnotationDirectoryFormat:
    return _collate(annotations)


def collate_mags_kmer_analyses(
    kmer_analyses: CARDMAGsKmerAnalysisDirectoryFormat,
) -> CARDMAGsKmerAnalysisDirectoryFormat:
    return _collate(kmer_analyses)


def collate_reads_allele_kmer_analyses(
    kmer_analyses: CARDReadsAlleleKmerAnalysisDirectoryFormat,
) -> CARDReadsAlleleKmerAnalysisDirectoryFormat:
    return _collate(kmer_analyses)


def collate_reads_gene_kmer_analyses(
    kmer_analyses: CARDReadsGeneKmerAnalysisDirectoryFormat,
) -> CARDReadsGeneKmerAnalysisDirectoryFormat:
    return _collate(kmer_analyses)


def _collate(partition_list):
    collated_partitions = type(partition_list[0])()
    # For every partition
    for partition in partition_list:
        # For every sample
        for sample in partition.path.iterdir():
            # If artifacts are annotations or kmer analyses from MAGs
            if isinstance(
                partition_list[0],
                (CARDAnnotationDirectoryFormat, CARDMAGsKmerAnalysisDirectoryFormat),
            ):
                # For every MAG
                for mag in sample.iterdir():
                    # Create directories in collate. If dir already exists raise error
                    try:
                        os.makedirs(collated_partitions.path / sample.name / mag.name)
                    except FileExistsError as e:
                        raise FileExistsError(
                            f"The directory already exists: {e.filename}. MAG IDs must"
                            f" be unique across all artifacts. Each artifact in the"
                            f" list must be unique and cannot be repeated."
                        )

                    # Copy every file in the MAG directory to the collated directory
                    for file in mag.iterdir():
                        duplicate(
                            file,
                            collated_partitions.path
                            / sample.name
                            / mag.name
                            / file.name,
                        )

            # If artifacts are annotations or kmer analyses are from reads
            else:
                # Create directories in collate. If dir already exists raise error
                try:
                    os.makedirs(collated_partitions.path / sample.name)
                except FileExistsError as e:
                    raise FileExistsError(
                        f"The directory already exists: {e.filename}. Sample IDs must"
                        f" be unique across all artifacts. Each artifact in the"
                        f" list must be unique and cannot be repeated."
                    )

                # Copy every file in the sample directory to the collated directory
                for file in sample.iterdir():
                    duplicate(file, collated_partitions.path / sample.name / file.name)

    return collated_partitions


def partition_mags_annotations(
    annotations: CARDAnnotationDirectoryFormat, num_partitions: int = None
) -> CARDAnnotationDirectoryFormat:
    return _partition_annotations(annotations, num_partitions)


def partition_reads_allele_annotations(
    annotations: CARDAlleleAnnotationDirectoryFormat, num_partitions: int = None
) -> CARDAlleleAnnotationDirectoryFormat:
    return _partition_annotations(annotations, num_partitions)


def partition_reads_gene_annotations(
    annotations: CARDGeneAnnotationDirectoryFormat, num_partitions: int = None
) -> CARDGeneAnnotationDirectoryFormat:
    return _partition_annotations(annotations, num_partitions)


def _partition_annotations(
    annotations: Union[
        CARDAnnotationDirectoryFormat,
        CARDGeneAnnotationDirectoryFormat,
        CARDAlleleAnnotationDirectoryFormat,
    ],
    num_partitions: int = None,
):
    partitioned_annotations = {}

    # Get dict with paths to all files in artifact and get number of annotations
    annotations_dict = annotations.sample_dict()
    num_annotations = len(annotations_dict)

    # If no number of partitions is specified or the number is higher than the number
    # of annotations, all annotations get partitioned by annotation
    if num_partitions is None:
        num_partitions = num_annotations
    elif num_partitions > num_annotations:
        warnings.warn(
            "You have requested a number of partitions"
            f" '{num_partitions}' that is greater than your number"
            f" of annotations '{num_annotations}'. Your data will be"
            f" partitioned by annotation into '{num_annotations}'"
            " partitions."
        )
        num_partitions = num_annotations

    # Split the dict into a list of specified number of dicts
    i = itertools.cycle(range(num_partitions))
    annotations_split_list = [{} for _ in range(num_partitions)]
    for key, value in annotations_dict.items():
        annotations_split_list[next(i)][key] = value

    for i, annotations_split in enumerate(annotations_split_list, 1):
        # Creates directory with same format as input
        partitioned_annotation = type(annotations)()

        # Constructs paths to all annotation files and move them to the new partition
        # directories
        for sample_mag_id, file_path_list in annotations_split.items():
            for file_path in file_path_list:
                file_path_des = os.path.join(
                    partitioned_annotation.path,
                    sample_mag_id,
                    os.path.basename(file_path),
                )
                os.makedirs(os.path.dirname(file_path_des), exist_ok=True)
                duplicate(file_path, file_path_des)

        if isinstance(annotations, CARDAnnotationDirectoryFormat):
            partitioned_annotation_key = sample_mag_id.replace("/", "_")
        else:
            partitioned_annotation_key = sample_mag_id

        # Add the partitioned object to the collection
        if num_partitions == num_annotations:  # and not duplicates:
            partitioned_annotations[partitioned_annotation_key] = partitioned_annotation
        else:
            partitioned_annotations[i] = partitioned_annotation

    return partitioned_annotations
