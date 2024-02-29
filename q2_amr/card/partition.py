import os

from qiime2.util import duplicate

from q2_amr.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
)


def collate_mags_annotations(
    annotations: CARDAnnotationDirectoryFormat,
) -> CARDAnnotationDirectoryFormat:
    return _collate_annotations(annotations)


def collate_reads_allele_annotations(
    annotations: CARDAlleleAnnotationDirectoryFormat,
) -> CARDAlleleAnnotationDirectoryFormat:
    return _collate_annotations(annotations)


def collate_reads_gene_annotations(
    annotations: CARDGeneAnnotationDirectoryFormat,
) -> CARDGeneAnnotationDirectoryFormat:
    return _collate_annotations(annotations)


def _collate_annotations(annotations):
    collated_annotations = type(annotations)()
    # For every partition
    for annotation in annotations:
        # For every sample
        for sample in annotation.path.iterdir():
            # If annotations are from MAGs
            if type(annotations) is CARDAnnotationDirectoryFormat:
                # For every MAG
                for mag in sample.iterdir():
                    # Create directories in collate
                    os.makedirs(
                        collated_annotations.path / sample.name / mag.name,
                        exist_ok=True,
                    )

                    # Copy every file in the MAG directory to the collated directory
                    for file in mag.iterdir():
                        duplicate(
                            file,
                            collated_annotations.path
                            / sample.name
                            / mag.name
                            / file.name,
                        )

            # If annotations are from reads
            else:
                # Create directories in collate
                os.makedirs(collated_annotations.path / sample.name, exist_ok=True)

                # For every mag in the sample
                for file in sample.iterdir():
                    duplicate(file, collated_annotations.path / sample.name / file.name)

    return collated_annotations
