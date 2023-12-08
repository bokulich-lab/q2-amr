import json
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


def load_preprocess_card_db(tmp, card_db, operation):
    if operation == "load":
        cmd = ["rgi", "load", "--card_json", str(card_db), "--local"]
    elif operation == "preprocess":
        cmd = ["rgi", "card_annotation", "-i", str(card_db)]
    elif operation == "load_fasta":
        with open(str(card_db)) as f:
            card_data = json.load(f)
            version = card_data["_version"]
        cmd = [
            "rgi",
            "load",
            "-i",
            str(card_db),
            "--card_annotation",
            f"card_database_v{version}.fasta",
            "--local",
        ]
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
