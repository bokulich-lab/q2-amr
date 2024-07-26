import os
import shutil
import subprocess
import tempfile

import pandas as pd
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt

from q2_amr.card.types import CARDAnnotationDirectoryFormat, CARDDatabaseDirectoryFormat
from q2_amr.card.utils import load_card_db
from q2_amr.utils.utils import create_count_table, read_in_txt, run_command


def annotate_mags_card(
    mag: MultiMAGSequencesDirFmt,
    card_db: CARDDatabaseDirectoryFormat,
    alignment_tool: str = "BLAST",
    split_prodigal_jobs: bool = False,
    include_loose: bool = False,
    include_nudge: bool = False,
    low_quality: bool = False,
    threads: int = 1,
) -> (CARDAnnotationDirectoryFormat, pd.DataFrame):
    manifest = mag.manifest.view(pd.DataFrame)
    amr_annotations = CARDAnnotationDirectoryFormat()
    frequency_list = []
    with tempfile.TemporaryDirectory() as tmp:
        load_card_db(card_db=card_db)
        for samp_bin in list(manifest.index):
            bin_dir = os.path.join(str(amr_annotations), samp_bin[0], samp_bin[1])
            os.makedirs(bin_dir, exist_ok=True)
            input_sequence = manifest.loc[samp_bin, "filename"]
            run_rgi_main(
                tmp,
                input_sequence,
                alignment_tool,
                split_prodigal_jobs,
                include_loose,
                include_nudge,
                low_quality,
                threads,
            )
            txt_path = os.path.join(bin_dir, "amr_annotation.txt")
            json_path = os.path.join(bin_dir, "amr_annotation.json")

            shutil.move(f"{tmp}/output.txt", txt_path)
            shutil.move(f"{tmp}/output.json", json_path)
            samp_bin_name = os.path.join(samp_bin[0], samp_bin[1])
            frequency_df = read_in_txt(
                path=txt_path,
                samp_bin_name=samp_bin_name,
                data_type="mags",
                colname="Best_Hit_ARO",
            )
            frequency_list.append(frequency_df)
        feature_table = create_count_table(df_list=frequency_list)
    return (
        amr_annotations,
        feature_table,
    )


def run_rgi_main(
    tmp,
    input_sequence: str,
    alignment_tool: str = "BLAST",
    split_prodigal_jobs: bool = False,
    include_loose: bool = False,
    include_nudge: bool = False,
    low_quality: bool = False,
    num_threads: int = 1,
):
    cmd = [
        "rgi",
        "main",
        "--input_sequence",
        input_sequence,
        "--output_file",
        f"{tmp}/output",
        "-n",
        str(num_threads),
        "--alignment_tool",
        alignment_tool,
        "--input_type",
        "contig",
    ]
    if include_loose:
        cmd.append("--include_loose")
    if include_nudge:
        cmd.append("--include_nudge")
    if low_quality:
        cmd.append("--low_quality")
    if split_prodigal_jobs:
        cmd.append("--split_prodigal_jobs")
    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
