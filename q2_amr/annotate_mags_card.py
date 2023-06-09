import os
import shutil
import subprocess
import tempfile

import pandas as pd
from q2_types_genomics.per_sample_data import MultiMAGSequencesDirFmt

from q2_amr.types import CARDAnnotationDirectoryFormat, CARDDatabaseFormat
from q2_amr.utils import load_preprocess_card_db, run_command


def annotate_mags_card(
    mag: MultiMAGSequencesDirFmt,
    card_db: CARDDatabaseFormat,
    alignment_tool: str = "BLAST",
    input_type: str = "contig",
    split_prodigal_jobs: bool = False,
    include_loose: bool = False,
    include_nudge: bool = False,
    low_quality: bool = False,
    num_threads: int = 1,
) -> CARDAnnotationDirectoryFormat:
    manifest = mag.manifest.view(pd.DataFrame)
    amr_annotations = CARDAnnotationDirectoryFormat()
    with tempfile.TemporaryDirectory() as tmp:
        load_preprocess_card_db(tmp, card_db, "load")
        for samp_bin in list(manifest.index):
            bin_dir = os.path.join(str(amr_annotations), samp_bin[0], samp_bin[1])
            os.makedirs(bin_dir, exist_ok=True)
            input_sequence = manifest.loc[samp_bin, "filename"]
            run_rgi_main(
                tmp,
                input_sequence,
                alignment_tool,
                input_type,
                split_prodigal_jobs,
                include_loose,
                include_nudge,
                low_quality,
                num_threads,
            )
            shutil.move(f"{tmp}/output.txt", f"{bin_dir}/amr_annotation.txt")
            shutil.move(f"{tmp}/output.json", f"{bin_dir}/amr_annotation.json")
    print("a")
    return amr_annotations


def run_rgi_main(
    tmp,
    input_sequence: str,
    alignment_tool: str = "BLAST",
    input_type: str = "contig",
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
        input_type,
        "--local",
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
