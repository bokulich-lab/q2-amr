import os
import shutil
import warnings
from typing import Union

import numpy as np

from q2_amr.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
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
    # Save all dir paths and file names of the annotations as dictionaries in a list
    annotations_all = []

    for dirpath, dirnames, filenames in os.walk(annotations.path):
        # This makes sure the location is the directory with the annotation files
        if not dirnames:
            components = os.path.normpath(dirpath).split(os.path.sep)
            dirs_and_files = {
                "dir_path": dirpath,
                "path_component_1": components[-1],
                "path_component_2": components[-2],
                "files": filenames,
            }
            annotations_all.append(dirs_and_files)

    # Retrieve the number of MAGs or reads annotations
    num_annotations = len(annotations_all)

    # If no number of partitions is specified or the number is higher than the number
    # of annotations, all annotations get partitioned into the number of annotations
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
    partition_dict = np.array_split(annotations_all, num_partitions)

    # Check if there are duplicates in the sample or MAG ids
    sample_mag_ids = [entry["path_component_1"] for entry in annotations_all]
    duplicates = True if len(sample_mag_ids) != len(set(sample_mag_ids)) else False

    for i, _partition_dict in enumerate(partition_dict, 1):
        # Creates directory with same format as input
        partitioned_annotation = type(annotations)()

        # Constructs paths to all annotation files and move them to the new partitioned
        # directories
        for dirs_and_files in _partition_dict:
            if type(annotations) is CARDAnnotationDirectoryFormat:
                result_dir = os.path.join(
                    partitioned_annotation.path,
                    dirs_and_files["path_component_2"],
                    dirs_and_files["path_component_1"],
                )
            else:
                result_dir = os.path.join(
                    partitioned_annotation.path, dirs_and_files["path_component_1"]
                )
            os.makedirs(result_dir)

            for file in dirs_and_files["files"]:
                shutil.copy(
                    os.path.join(dirs_and_files["dir_path"], file),
                    os.path.join(result_dir, file),
                )

            # Adds the partitioned object to the collection
            if num_partitions == num_annotations and not duplicates:
                partitioned_annotations[
                    dirs_and_files["path_component_1"]
                ] = partitioned_annotation
            else:
                partitioned_annotations[i] = partitioned_annotation

    return partitioned_annotations
