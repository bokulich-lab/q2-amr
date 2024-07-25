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


def run_command(cmd, cwd=None, verbose=True):
    if verbose:
        print(EXTERNAL_CMD_WARNING)
        print("\nCommand:", end=" ")
        print(" ".join(cmd), end="\n\n")
    subprocess.run(cmd, check=True, cwd=cwd)


def colorify(string: str):
    return "%s%s%s" % ("\033[1;32m", string, "\033[0m")


def read_in_txt(path: str, samp_bin_name: str, data_type: str, colname: str):
    # Read in txt file to pd.Dataframe
    df = pd.read_csv(path, sep="\t")

    # Process the df depending on the data type
    if data_type == "reads":
        df = df[[colname, "All Mapped Reads"]]
        df.rename(columns={"All Mapped Reads": samp_bin_name}, inplace=True)
    else:
        df = df[colname].value_counts().reset_index()
        df.columns = [colname, samp_bin_name]

    df = df.astype(str)
    return df


def create_count_table(df_list: list) -> pd.DataFrame:
    # Remove all empty lists from df_list
    df_list = [df for df in df_list if not df.empty]

    # Raise ValueError if df_list is empty. This happens when no ARGs were detected
    if not df_list:
        raise ValueError("No AMR genes where identified. No output can be created.")

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


def create_append_df(file_path, df_list, id_name, id_value):
    # Read in df
    df = pd.read_csv(file_path, sep="\t")

    # Insert column with sample or mag IDs
    df.insert(0, id_name, id_value)

    # Append df to df list
    df_list.append(df)


def combine_dataframes(df_list):
    # Concat all dfs
    df_combined = pd.concat(df_list, axis=0)

    # Sort all values by sample/mag ID column
    df_combined.sort_values(by=df_combined.columns[0], inplace=True)

    # Reset and rename index and set it to string to conform to metadata format
    df_combined.reset_index(inplace=True, drop=True)
    df_combined.index.name = "id"
    df_combined.index = df_combined.index.astype(str)

    return df_combined
