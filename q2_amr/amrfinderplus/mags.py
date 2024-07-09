import os
import shutil
import tempfile

import pandas as pd
from q2_types.genome_data import GenesDirectoryFormat
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt

from q2_amr.amrfinderplus.types import (
    AMRFinderPlusDatabaseDirFmt,
    ARMFinderPlusAnnotationsDirFmt,
)
from q2_amr.amrfinderplus.utils import run_amrfinderplus_n
from q2_amr.card.utils import create_count_table, read_in_txt


def annotate_mags_amrfinderplus(
    mags: MultiMAGSequencesDirFmt,
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
    manifest = mags.manifest.view(pd.DataFrame)

    annotations = ARMFinderPlusAnnotationsDirFmt()
    mutations = ARMFinderPlusAnnotationsDirFmt()
    genes = GenesDirectoryFormat()

    frequency_list = []

    with tempfile.TemporaryDirectory() as tmp:
        for samp_mag in list(manifest.index):
            input_sequence = manifest.loc[samp_mag, "filename"]

            sample_id = samp_mag[0]
            mag_id = samp_mag[1]

            run_amrfinderplus_n(
                working_dir=tmp,
                amrfinderplus_db=amrfinderplus_db,
                dna_sequence=input_sequence,
                protein_sequence=None,
                gff=None,
                organism=organism,
                plus=plus,
                report_all_equal=report_all_equal,
                ident_min=ident_min,
                coverage_min=coverage_min,
                translation_table=translation_table,
                threads=threads,
                id=mag_id + "_",
            )

            frequency_df = read_in_txt(
                path=os.path.join(tmp, f"{mag_id}_amr_annotations.tsv"),
                samp_bin_name=str(os.path.join(sample_id, mag_id)),
                data_type="mags",
                colname="Gene symbol",
            )

            for dir_format, file_name in zip(
                [annotations, mutations, genes],
                [
                    f"{mag_id}_amr_annotations.tsv",
                    f"{mag_id}_amr_mutations.tsv",
                    f"{mag_id}_amr_genes.fasta",
                ],
            ):
                if dir_format in [annotations, mutations]:
                    des_dir = os.path.join(str(dir_format), sample_id)
                    os.makedirs(des_dir, exist_ok=True)
                else:
                    des_dir = str(dir_format)
                shutil.move(os.path.join(tmp, file_name), des_dir)

            frequency_list.append(frequency_df)

        feature_table = create_count_table(df_list=frequency_list)
    return (
        annotations,
        mutations,
        genes,
        feature_table,
    )
