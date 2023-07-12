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
