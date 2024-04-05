import os
import shutil
import subprocess
import tempfile
from typing import Union

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)

from q2_amr.card.utils import create_count_table, load_card_db, read_in_txt, run_command
from q2_amr.types import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDDatabaseDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
)


def annotate_reads_card(
    reads: Union[
        SingleLanePerSamplePairedEndFastqDirFmt, SingleLanePerSampleSingleEndFastqDirFmt
    ],
    card_db: CARDDatabaseDirectoryFormat,
    aligner: str = "kma",
    threads: int = 1,
    include_wildcard: bool = False,
    include_other_models: bool = False,
) -> (
    CARDAlleleAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
    pd.DataFrame,
    pd.DataFrame,
):
    paired = isinstance(reads, SingleLanePerSamplePairedEndFastqDirFmt)
    manifest = reads.manifest.view(pd.DataFrame)
    allele_frequency_list, gene_frequency_list = [], []
    amr_allele_annotation = CARDAlleleAnnotationDirectoryFormat()
    amr_gene_annotation = CARDGeneAnnotationDirectoryFormat()
    with tempfile.TemporaryDirectory() as tmp:
        load_card_db(
            tmp=tmp,
            card_db=card_db,
            fasta=True,
            include_other_models=include_other_models,
            include_wildcard=include_wildcard,
        )
        for samp in list(manifest.index):
            fwd = manifest.loc[samp, "forward"]
            rev = manifest.loc[samp, "reverse"] if paired else None
            samp_allele_dir = os.path.join(str(amr_allele_annotation), samp)
            samp_gene_dir = os.path.join(str(amr_gene_annotation), samp)
            os.makedirs(samp_allele_dir)
            os.makedirs(samp_gene_dir)
            samp_input_dir = os.path.join(tmp, samp)
            os.makedirs(samp_input_dir)
            run_rgi_bwt(
                cwd=tmp,
                samp=samp,
                fwd=fwd,
                rev=rev,
                aligner=aligner,
                threads=threads,
                include_wildcard=include_wildcard,
                include_other_models=include_other_models,
            )

            # Create a frequency table and add it to a list, for gene and allele
            # mapping data
            for map_type, table_list in zip(
                ["allele", "gene"], [allele_frequency_list, gene_frequency_list]
            ):
                path_txt = os.path.join(
                    samp_input_dir, f"output.{map_type}_mapping_data.txt"
                )
                frequency_table = read_in_txt(
                    path=path_txt,
                    samp_bin_name=samp,
                    data_type="reads",
                    map_type=map_type,
                )
                table_list.append(frequency_table)

            # Move mapping and stats files to the sample allele and gene directories
            for map_type, des_dir in zip(
                ["allele", "gene"], [samp_allele_dir, samp_gene_dir]
            ):
                files = [f"{map_type}_mapping_data.txt"]
                # mapping statistics only go to the allele directories
                files.extend(
                    ["overall_mapping_stats.txt", "sorted.length_100.bam"]
                ) if map_type == "allele" else None
                for file in files:
                    shutil.copy(
                        os.path.join(samp_input_dir, "output." + file),
                        os.path.join(des_dir, file),
                    )

    allele_feature_table = create_count_table(allele_frequency_list)
    gene_feature_table = create_count_table(gene_frequency_list)
    return (
        amr_allele_annotation,
        amr_gene_annotation,
        allele_feature_table,
        gene_feature_table,
    )


def run_rgi_bwt(
    cwd: str,
    samp: str,
    fwd: str,
    rev: str,
    aligner: str,
    threads: int,
    include_wildcard: bool,
    include_other_models: bool,
):
    cmd = [
        "rgi",
        "bwt",
        "--read_one",
        fwd,
        "--output_file",
        f"{cwd}/{samp}/output",
        "-n",
        str(threads),
        "--clean",
        "--aligner",
        aligner,
    ]
    if rev:
        cmd.extend(["--read_two", rev])
    if include_wildcard:
        cmd.append("--include_wildcard")
    if include_other_models:
        cmd.append("--include_other_models")
    try:
        run_command(cmd, cwd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
