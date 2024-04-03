import os
import warnings
from typing import Union

import numpy as np
from qiime2.util import duplicate

from q2_amr.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
)


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
    annotations_all = []

    # Add one tuple with sample id MAG id and full paths to annotation files to
    # annotations_all per annotation file
    if isinstance(annotations, CARDAnnotationDirectoryFormat):
        for sample_id, mag in annotations.sample_dict().items():
            for mag_id, annotation_fp_list in mag.items():
                for annotation_fp in annotation_fp_list:
                    annotations_all.append((sample_id, mag_id, annotation_fp))

    else:
        for sample_id, annotation_fp_list in annotations.sample_dict().items():
            for annotation_fp in annotation_fp_list:
                annotations_all.append((sample_id, annotation_fp))

    # Retrieve the number of annotations
    num_annotations = len({tup[-2] for tup in annotations_all})

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

    # Splits annotations into the specified number of partitions
    arrays = np.array_split(annotations_all, num_partitions)

    for i, annotation_tuple in enumerate(arrays, 1):
        # Creates directory with same format as input
        partitioned_annotation = type(annotations)()

        # Constructs paths to all annotation files and move them to the new partition
        # directories
        if isinstance(annotations, CARDAnnotationDirectoryFormat):
            for sample_id, mag_id, annotation_fp in annotation_tuple:
                annotation_des_fp = os.path.join(
                    partitioned_annotation.path,
                    sample_id,
                    mag_id,
                    os.path.basename(annotation_fp),
                )
                os.makedirs(os.path.dirname(annotation_des_fp), exist_ok=True)
                duplicate(annotation_fp, annotation_des_fp)

                partitioned_annotation_key = mag_id

        else:
            for sample_id, annotation_fp in annotation_tuple:
                annotation_des_fp = os.path.join(
                    partitioned_annotation.path,
                    sample_id,
                    os.path.basename(annotation_fp),
                )
                os.makedirs(os.path.dirname(annotation_des_fp), exist_ok=True)
                duplicate(annotation_fp, annotation_des_fp)

                partitioned_annotation_key = sample_id

        # Add the partitioned object to the collection
        if num_partitions == num_annotations:  # and not duplicates:
            partitioned_annotations[partitioned_annotation_key] = partitioned_annotation
        else:
            partitioned_annotations[i] = partitioned_annotation

    return partitioned_annotations


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
    for annotation in partition_list:
        # For every sample
        for sample in annotation.path.iterdir():
            # If formats are annotations or kmer analyses from MAGs
            if isinstance(
                partition_list[0],
                (CARDAnnotationDirectoryFormat, CARDMAGsKmerAnalysisDirectoryFormat),
            ):
                # For every MAG
                for mag in sample.iterdir():
                    # Create directories in collate
                    os.makedirs(
                        collated_partitions.path / sample.name / mag.name,
                        exist_ok=True,
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

            # If annotations or kmer analyses are from reads
            else:
                # Create directories in collate object
                os.makedirs(collated_partitions.path / sample.name, exist_ok=True)

                # For every mag in the sample
                for file in sample.iterdir():
                    duplicate(file, collated_partitions.path / sample.name / file.name)

    return collated_partitions
