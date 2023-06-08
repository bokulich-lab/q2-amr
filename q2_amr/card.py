import os
import shutil
import subprocess
import tarfile
import tempfile
from distutils.dir_util import copy_tree
from functools import reduce
from typing import Union

import altair as alt
import pandas as pd
import pkg_resources
import q2templates
import requests
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from q2_types_genomics.per_sample_data import MultiMAGSequencesDirFmt

from q2_amr.types import (
    CARDAnnotationDirectoryFormat,
    CARDAnnotationJSONFormat,
    CARDDatabaseDirectoryFormat,
    CARDDatabaseFormat,
)
from q2_amr.types._format import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
)
from q2_amr.utils import load_preprocess_card_db, run_command

CARD_URL = "https://card.mcmaster.ca/latest/data"


def fetch_card_db() -> CARDDatabaseDirectoryFormat:
    try:
        response = requests.get(CARD_URL, stream=True)
    except requests.ConnectionError as e:
        raise requests.ConnectionError("Network connectivity problems.") from e
    with tempfile.TemporaryDirectory() as tmp_dir:
        try:
            with tarfile.open(fileobj=response.raw, mode="r|bz2") as tar:
                tar.extractall(path=tmp_dir)
        except tarfile.ReadError as a:
            raise tarfile.ReadError("Tarfile is invalid.") from a
        card_db = CARDDatabaseDirectoryFormat()
        shutil.move(
            os.path.join(tmp_dir, "card.json"), os.path.join(str(card_db), "card.json")
        )
        return card_db


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


def heatmap(
    output_dir: str,
    amr_annotation_json: CARDAnnotationJSONFormat,
    # clus: str = 'no',
    # cat: str = 'no',
    # frequency=False
) -> None:
    TEMPLATES = pkg_resources.resource_filename("q2_amr", "assets")
    with tempfile.TemporaryDirectory() as tmp:
        results_dir = os.path.join(tmp, "results")
        os.makedirs(results_dir)
        cmd = [
            f"rgi heatmap --input {os.path.dirname(str(amr_annotation_json))} --output "
            f"{tmp}/results/heatmap"
        ]
        # if frequency:
        #     cmd.extend(["--frequency"])
        # if clus == 'both':
        #     cmd.extend(["-clus both"])
        # elif clus == 'samples':
        #     cmd.extend(["-clus samples"])
        # elif clus == 'genes':
        #     cmd.extend(["-clus genes"])
        # elif clus == 'no':
        #     pass
        # if cat == 'drug_class':
        #     cmd.extend(["-cat drug_class"])
        # elif cat == 'resistance_mechanism':
        #     cmd.extend(["-cat resistance_mechanism"])
        # elif cat == 'gene_family':
        #     cmd.extend(["-cat gene_family"])
        # elif cat == 'no':
        #     pass
        try:
            run_command(cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running rgi, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )
        extensions = [
            ".eps",
            ".csv",
            ".png",
        ]  # Replace with the file extensions you want to rename
        # Get all the files in the directory
        files = os.listdir(results_dir)
        # Loop through all the files in the directory
        for filename in files:
            if os.path.splitext(filename)[1] in extensions:
                # Construct the new file name with the correct extension
                file_ext = os.path.splitext(filename)[1]
                new_filename = "heatmap" + file_ext

                # Rename the file
                old_path = os.path.join(results_dir, filename)
                new_path = os.path.join(results_dir, new_filename)
                os.rename(old_path, new_path)
        copy_tree(os.path.join(TEMPLATES, "rgi"), output_dir)
        copy_tree(results_dir, os.path.join(output_dir, "rgi_data"))
    context = {"tabs": [{"title": "QC report", "url": "index.html"}]}
    index = os.path.join(TEMPLATES, "rgi", "index.html")
    templates = [index]
    q2templates.render(templates, output_dir, context=context)


def annotate_reads_card(
    reads: Union[
        SingleLanePerSamplePairedEndFastqDirFmt, SingleLanePerSampleSingleEndFastqDirFmt
    ],
    card_db: CARDDatabaseFormat,
    aligner: str = "kma",
    threads: int = 1,
    include_baits: bool = False,
    mapq: float = None,
    mapped: float = None,
    coverage: float = None,
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
        load_preprocess_card_db(tmp, card_db, "load")
        load_preprocess_card_db(tmp, card_db, "preprocess")
        load_preprocess_card_db(tmp, card_db, "load_fasta")
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
                include_baits=include_baits,
                mapq=mapq,
                mapped=mapped,
                coverage=coverage,
            )
            allele_frequency = read_in_txt(samp_input_dir, "allele")
            allele_frequency_list.append(allele_frequency)
            gene_frequency = read_in_txt(samp_input_dir, "gene")
            gene_frequency_list.append(gene_frequency)
            move_files(samp_input_dir, samp_allele_dir, "allele")
            move_files(samp_input_dir, samp_gene_dir, "gene")

    allele_feature_table = create_count_table(allele_frequency_list)
    gene_feature_table = create_count_table(gene_frequency_list)
    return (
        amr_allele_annotation,
        amr_gene_annotation,
        allele_feature_table,
        gene_feature_table,
    )


def move_files(source_dir: str, des_dir: str, map_type: str):
    shutil.move(
        os.path.join(source_dir, f"output.{map_type}_mapping_data.txt"),
        os.path.join(des_dir, f"{map_type}_mapping_data.txt"),
    )
    shutil.copy(
        os.path.join(source_dir, "output.overall_mapping_stats.txt"),
        os.path.join(des_dir, "overall_mapping_stats.txt"),
    )


def create_count_table(df_list: list) -> pd.DataFrame:
    df_merged = reduce(
        lambda left, right: pd.merge(left, right, on="ARO Accession", how="outer"),
        df_list,
    )
    df_transposed = df_merged.transpose()
    df_transposed = df_transposed.fillna(0)
    df_transposed.columns = df_transposed.iloc[0]
    df_transposed = df_transposed.drop("ARO Accession")
    df_transposed.columns.name = None
    df_transposed.index.name = "sample_id"
    return df_transposed


def read_in_txt(samp_dir: str, map_type: str):
    df = pd.read_csv(
        os.path.join(samp_dir, f"output.{map_type}_mapping_data.txt"), sep="\t"
    )
    df = df[["ARO Accession"]]
    df = df.astype(str)
    samp = os.path.basename(samp_dir)
    df[samp] = df.groupby("ARO Accession")["ARO Accession"].transform("count")
    df = df.drop_duplicates(subset=["ARO Accession"])
    return df


def run_rgi_bwt(
    cwd: str,
    samp: str,
    fwd: str,
    rev: str,
    aligner: str,
    threads: int,
    include_baits: bool,
    mapq: float,
    mapped: float,
    coverage: float,
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
    if include_baits:
        cmd.append("--include_baits")
    if mapq:
        cmd.extend(["--mapq", str(mapq)])
    if mapped:
        cmd.extend(["--mapped", str(mapped)])
    if coverage:
        cmd.extend(["--coverage", str(coverage)])
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
