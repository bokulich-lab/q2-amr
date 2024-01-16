import os
import shutil
import subprocess
import tempfile
from distutils.dir_util import copy_tree
from typing import Union

import altair as alt
import pandas as pd
import pkg_resources
import q2templates
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
                files.append(
                    "overall_mapping_stats.txt"
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
        "--local",
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


def plot_sample_stats(sample_stats: dict, output_dir: str):
    sample_stats_df = pd.DataFrame.from_dict(sample_stats, orient="index")
    sample_stats_df.reset_index(inplace=True)

    mapped_reads_plot = (
        alt.Chart(sample_stats_df)
        .mark_bar()
        .encode(
            x=alt.X(
                "index:N",
                title=None,
                axis=alt.Axis(labels=False, ticks=False),
                bandPosition=0.1,
            ),
            y=alt.Y(
                "mapped_reads",
                title="Mapped Reads",
                scale=alt.Scale(
                    domain=(0, sample_stats_df["mapped_reads"].max() * 1.1)
                ),
            ),
            tooltip=[
                alt.Tooltip("index", title="Sample"),
                alt.Tooltip("mapped_reads", title="Mapped Reads"),
                alt.Tooltip("total_reads", title="Total Reads"),
            ],
        )
        .properties(width=alt.Step(80), height=200)
    )

    percentage_plot = (
        alt.Chart(sample_stats_df)
        .mark_bar()
        .encode(
            x=alt.X(
                "index:N", title=None, axis=alt.Axis(labelAngle=15), bandPosition=0.1
            ),
            y=alt.Y(
                "percentage",
                title="Mapped Reads (%)",
                scale=alt.Scale(domain=(0, sample_stats_df["percentage"].max() * 1.1)),
            ),
            tooltip=[
                alt.Tooltip("index", title="Sample"),
                alt.Tooltip("percentage", title="Mapped Reads (%)"),
                alt.Tooltip("total_reads", title="Total Reads"),
            ],
        )
        .properties(width=alt.Step(80), height=200)
    )
    combined_chart = alt.vconcat(mapped_reads_plot, percentage_plot, spacing=0)
    combined_chart.save(os.path.join(output_dir, "sample_stats_plot.html"))


def extract_sample_stats(samp_dir: str):
    with open(os.path.join(samp_dir, "overall_mapping_stats.txt"), "r") as f:
        for line in f:
            if "Total reads:" in line:
                total_reads = int(line.split()[2])
            elif "Mapped reads:" in line:
                mapped_reads = int(line.split()[2])
                percentage = float(line.split()[3].strip("()").strip("%"))
        sample_stats_dict = {
            "total_reads": total_reads,
            "mapped_reads": mapped_reads,
            "percentage": percentage,
        }
    return sample_stats_dict


def visualize_annotation_stats(
    output_dir: str,
    amr_reads_annotation: Union[
        CARDGeneAnnotationDirectoryFormat, CARDAlleleAnnotationDirectoryFormat
    ],
):
    directory = str(amr_reads_annotation)
    sample_stats = {}
    for samp in os.listdir(directory):
        samp_dir = os.path.join(directory, samp)
        sample_stats[samp] = extract_sample_stats(samp_dir)
    plot_sample_stats(sample_stats, output_dir)
    TEMPLATES = pkg_resources.resource_filename("q2_amr", "assets")
    copy_tree(os.path.join(TEMPLATES, "rgi", "annotation_stats"), output_dir)
    context = {"tabs": [{"title": "Mapped Reads", "url": "index.html"}]}
    index = os.path.join(TEMPLATES, "rgi", "annotation_stats", "index.html")
    templates = [index]
    q2templates.render(templates, output_dir, context=context)
