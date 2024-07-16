import os
import shutil
import tempfile
from typing import Union

import pandas as pd
from q2_types.genome_data import GenesDirectoryFormat
from q2_types.per_sample_sequences import ContigSequencesDirFmt, MultiMAGSequencesDirFmt

from q2_amr.amrfinderplus.types import (
    AMRFinderPlusAnnotationsDirFmt,
    AMRFinderPlusDatabaseDirFmt,
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
    curated_ident: bool = False,
    coverage_min: float = 0.5,
    translation_table: str = "11",
    threads: int = None,
) -> (
    AMRFinderPlusAnnotationsDirFmt,
    AMRFinderPlusAnnotationsDirFmt,
    GenesDirectoryFormat,
    pd.DataFrame,
):
    annotations = AMRFinderPlusAnnotationsDirFmt()
    mutations = AMRFinderPlusAnnotationsDirFmt()
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
            else:
                sample_id = os.path.splitext(os.path.basename(file))[0][:-8]
                mag_id = ""

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
                curated_ident=curated_ident,
                coverage_min=coverage_min,
                translation_table=translation_table,
                threads=threads,
            )

            # Create frequency dataframe and append it to list
            frequency_df = read_in_txt(
                path=os.path.join(tmp, "amr_annotations.tsv"),
                samp_bin_name=str(os.path.join(sample_id, mag_id)),
                data_type="mags",
                colname="Gene symbol",
            )
            frequency_list.append(frequency_df)

            # Move mutations file. If it is not created, create an empty mutations file
            des_path_mutations = os.path.join(
                str(mutations),
                sample_id,
                f"{mag_id + '_' if mag_id else ''}amr_mutations.tsv",
            )
            os.makedirs(os.path.dirname(des_path_mutations), exist_ok=True)
            if organism:
                shutil.move(os.path.join(tmp, "amr_mutations.tsv"), des_path_mutations)
            else:
                with open(des_path_mutations, "w"):
                    pass

            # Move annotations file
            des_path_annotations = os.path.join(
                str(annotations),
                sample_id,
                f"{mag_id + '_' if mag_id else ''}amr_annotations.tsv",
            )
            os.makedirs(os.path.dirname(des_path_annotations), exist_ok=True)
            shutil.move(os.path.join(tmp, "amr_annotations.tsv"), des_path_annotations)

            # Move genes file
            shutil.move(
                os.path.join(tmp, "amr_genes.fasta"),
                os.path.join(
                    str(genes), f"{mag_id if mag_id else sample_id}_amr_genes.fasta"
                ),
            )

        feature_table = create_count_table(df_list=frequency_list)
    return (
        annotations,
        mutations,
        genes,
        feature_table,
    )
