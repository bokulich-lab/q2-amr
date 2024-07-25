import pandas as pd


def create_df(file_path, df_list, id_name, id_value):
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
