import os
import shutil
import tempfile
from typing import Union

import pandas as pd
from q2_types.genome_data import GenesDirectoryFormat
from q2_types.per_sample_sequences import ContigSequencesDirFmt, MultiMAGSequencesDirFmt

from q2_amr.amrfinderplus.types import (
    AMRFinderPlusDatabaseDirFmt,
    ARMFinderPlusAnnotationsDirFmt,
)
from q2_amr.amrfinderplus.utils import run_amrfinderplus_n
from q2_amr.card.utils import create_count_table, read_in_txt


def annotate_sample_data_amrfinderplus(
    sequences: Union[MultiMAGSequencesDirFmt, ContigSequencesDirFmt],
    amrfinderplus_db: AMRFinderPlusDatabaseDirFmt,
    organism: str = None,
    plus: bool = False,
    report_all_equal: bool = False,
    ident_min: float = None,
    coverage_min: float = 0.5,
    translation_table: str = "11",
    threads: int = None,
) -> (
    ARMFinderPlusAnnotationsDirFmt,
    ARMFinderPlusAnnotationsDirFmt,
    GenesDirectoryFormat,
    pd.DataFrame,
):
    annotations = ARMFinderPlusAnnotationsDirFmt()
    mutations = ARMFinderPlusAnnotationsDirFmt()
    genes = GenesDirectoryFormat()
    frequency_list = []

    # Create list of paths to all mags or contigs
    if isinstance(sequences, MultiMAGSequencesDirFmt):
        manifest = sequences.manifest.view(pd.DataFrame)
        files = manifest["filename"]
    else:
        files = [
            os.path.join(str(sequences), file) for file in os.listdir(str(sequences))
        ]

    with tempfile.TemporaryDirectory() as tmp:
        # Iterate over paths of mags or contigs
        for file in files:
            # Set sample and mag ids and output file pats for mag or contig
            if isinstance(sequences, MultiMAGSequencesDirFmt):
                index_value = manifest.query("filename == @file").index[0]
                sample_id = index_value[0]
                mag_id = index_value[1]
                annotations_path = os.path.join(tmp, f"{mag_id}_amr_annotations.tsv")
                mutations_path = os.path.join(tmp, f"{mag_id}_amr_mutations.tsv")
                genes_path = os.path.join(tmp, f"{mag_id}_amr_genes.fasta")
            else:
                sample_id = os.path.splitext(os.path.basename(file))[0][:-8]
                mag_id = ""
                annotations_path = os.path.join(tmp, "amr_annotations.tsv")
                mutations_path = os.path.join(tmp, "amr_mutations.tsv")
                genes_path = os.path.join(tmp, f"{sample_id}_amr_genes.fasta")

            # Run amrfinderplus
            run_amrfinderplus_n(
                working_dir=tmp,
                amrfinderplus_db=amrfinderplus_db,
                dna_sequence=file,
                protein_sequence=None,
                gff=None,
                organism=organism,
                plus=plus,
                report_all_equal=report_all_equal,
                ident_min=ident_min,
                coverage_min=coverage_min,
                translation_table=translation_table,
                threads=threads,
                mag_id=mag_id,
                sample_id=sample_id,
            )

            # Create frequency dataframe and append it to list
            frequency_df = read_in_txt(
                path=os.path.join(tmp, annotations_path),
                samp_bin_name=str(os.path.join(sample_id, mag_id)),
                data_type="mags",
                colname="Gene symbol",
            )
            frequency_list.append(frequency_df)

            # Move mutations file. If it is not created, create an empty mutations file
            des_dir_mutations = os.path.join(str(mutations), sample_id)
            os.makedirs(des_dir_mutations, exist_ok=True)
            if organism:
                shutil.move(mutations_path, des_dir_mutations)
            else:
                with open(
                    os.path.join(
                        str(mutations),
                        des_dir_mutations,
                        os.path.basename(mutations_path),
                    ),
                    "w",
                ):
                    pass

            # Move annotations file
            des_dir_annotations = os.path.join(str(annotations), sample_id)
            os.makedirs(des_dir_annotations, exist_ok=True)
            shutil.move(annotations_path, des_dir_annotations)

            # Move genes file
            shutil.move(genes_path, str(genes))

        feature_table = create_count_table(df_list=frequency_list)
    return (
        annotations,
        mutations,
        genes,
        feature_table,
    )
