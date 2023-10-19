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


def load_card_db(tmp, card_db, operation, all_models, wildcard):
    path_card_json = os.path.join(str(card_db), "card.json")
    if operation == "load":
        cmd = ["rgi", "load", "--card_json", path_card_json, "--local"]
    elif operation == "load_fasta":
        with open(path_card_json) as f:
            card_data = json.load(f)
            version = card_data["_version"]
        models = ("_all", "_all_models") if all_models is True else ("", "")
        path_card_fasta = os.path.join(
            str(card_db), f"card_database_v{version}{models[0]}.fasta"
        )
        cmd = [
            "rgi",
            "load",
            "-i",
            path_card_json,
            f"--card_annotation{models[1]}",
            path_card_fasta,
            "--local",
        ]
        if wildcard:
            path_wildcard_fasta = os.path.join(
                str(card_db), f"wildcard_database_v0{models[0]}.fasta"
            )
            path_wildcard_index = os.path.join(
                str(card_db), "index-for-model-sequences.txt"
            )
            cmd.extend(
                [
                    f"--wildcard_annotation{models[1]}",
                    path_wildcard_fasta,
                    "--wildcard_index",
                    path_wildcard_index,
                ]
            )

    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def read_in_txt(path: str, col_name: str, samp_bin_name: str):
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return None
    df = df[[col_name]]
    df = df.astype(str)
    df[samp_bin_name] = df.groupby(col_name)[col_name].transform("count")
    df = df.drop_duplicates(subset=[col_name])
    return df


def create_count_table(df_list: list) -> pd.DataFrame:
    if not df_list:
        raise ValueError(
            "RGI did not identify any AMR genes. No output can be created."
        )
    df = reduce(
        lambda left, right: pd.merge(left, right, on=left.columns[0], how="outer"),
        df_list,
    )
    df = df.transpose()
    df = df.fillna(0)
    df.columns = df.iloc[0]
    df = df.drop(df.index[0])
    df.columns.name = None
    df.index.name = "sample_id"
    return df
