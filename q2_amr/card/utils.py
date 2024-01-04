import glob
import json
import os
import subprocess
from functools import reduce

import pandas as pd

EXTERNAL_CMD_WARNING = (
    "Running external command line application(s). "
    "This may print messages to stdout and/or stderr.\n"
    "The command(s) being run are below. These commands "
    "cannot be manually re-run as they will depend on "
    "temporary files that no longer exist."
)


def run_command(cmd, cwd, verbose=True):
    if verbose:
        print(EXTERNAL_CMD_WARNING)
        print("\nCommand:", end=" ")
        print(" ".join(cmd), end="\n\n")
    subprocess.run(cmd, check=True, cwd=cwd)


def load_card_db(
    tmp,
    card_db,
    kmer_db=None,
    kmer: bool = False,
    fasta: bool = False,
    include_other_models: bool = False,
    include_wildcard: bool = False,
):
    # Get path to card.json
    path_card_json = os.path.join(str(card_db), "card.json")

    # Base command that only loads card.json into the local database
    cmd = ["rgi", "load", "--card_json", path_card_json, "--local"]

    # Define suffixes for card fasta file
    models = ("_all", "_all_models") if include_other_models is True else ("", "")

    # Extend base command with flag to load card fasta file
    if fasta:
        # Retrieve the database version number from card.jason file
        with open(path_card_json) as f:
            card_data = json.load(f)
            version = card_data["_version"]

        # Define path to card fasta file
        path_card_fasta = os.path.join(
            str(card_db), f"card_database_v{version}{models[0]}.fasta"
        )

        # Extend base command
        cmd.extend([f"--card_annotation{models[1]}", path_card_fasta])

    # Extend base command with flag to load wildcard fasta file and index
    if include_wildcard:
        cmd.extend(
            [
                f"--wildcard_annotation{models[1]}",
                os.path.join(str(card_db), f"wildcard_database_v0{models[0]}.fasta"),
                "--wildcard_index",
                os.path.join(str(card_db), "index-for-model-sequences.txt"),
            ]
        )
    # Extend base command with flag to load kmer json and txt database files
    if kmer:
        path_kmer_json = glob.glob(os.path.join(str(kmer_db), "*_kmer_db.json"))[0]
        cmd.extend(
            [
                "--kmer_database",
                path_kmer_json,
                "--amr_kmers",
                glob.glob(os.path.join(str(kmer_db), "all_amr_*mers.txt"))[0],
                "--kmer_size",
                os.path.basename(path_kmer_json).split("_")[0],
            ]
        )

    # Run command
    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def read_in_txt(path: str, samp_bin_name: str, data_type):
    # Read in txt file to pd.Dataframe
    df = pd.read_csv(path, sep="\t")

    # Process the df depending on the data type (from reads or mags)
    if data_type == "reads":
        df = df[["ARO Term", "All Mapped Reads"]]
        df.rename(columns={"All Mapped Reads": samp_bin_name}, inplace=True)
    else:
        df = df[["Best_Hit_ARO"]]
        df[samp_bin_name] = 1

    df = df.astype(str)
    return df


def create_count_table(df_list: list) -> pd.DataFrame:
    # Remove all empty lists from df_list
    df_list = [df for df in df_list if not df.empty]

    # Raise ValueError if df_list is empty. This happens when no ARGs were detected
    if not df_list:
        raise ValueError(
            "RGI did not identify any AMR genes. No output can be created."
        )

    # Merge all dfs contained in df_list
    df = reduce(
        lambda left, right: pd.merge(left, right, on=left.columns[0], how="outer"),
        df_list,
    )

    # Process the df to meet all requirements for a FeatureTable
    df = df.transpose()
    df = df.fillna(0)
    df.columns = df.iloc[0]
    df = df.drop(df.index[0])
    df.columns.name = None
    df.index.name = "sample_id"
    return df
