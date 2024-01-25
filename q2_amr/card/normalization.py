import os
import subprocess
import tempfile

import biom
import pandas as pd
from pydeseq2.dds import DeseqDataSet

from q2_amr.card.utils import run_command
from q2_amr.types import CARDAlleleAnnotationDirectoryFormat


def normalize_mor(
    table: biom.Table,
) -> pd.DataFrame:
    # manifest = reads.manifest.view(pd.DataFrame)
    # #paired = isinstance(reads, SingleLanePerSamplePairedEndFastqDirFmt)
    # for samp in list(manifest.index):
    #     # fwd = manifest.loc[samp, "forward"]
    #     # rev = manifest.loc[samp, "reverse"] if paired else None
    data = table.matrix_data.toarray()
    columns = table.ids(axis="sample")
    index = table.ids(axis="observation")
    df = pd.DataFrame(data, index=index, columns=columns)
    df = df.T
    metadata = df.iloc[:, 0].to_frame()
    metadata.iloc[0, 0] = "a"
    original_columns = metadata.columns

    metadata.columns = ["condition"] + list(original_columns[1:])
    dds = DeseqDataSet(counts=df, metadata=metadata)
    dds.fit_size_factors()

    # normalizedcounts = dds.obsm["size_factors"]
    return df


def normalize_tpm(
    table: biom.Table, amr_annotations: CARDAlleleAnnotationDirectoryFormat
) -> pd.DataFrame:
    # Initialize an empty list to store individual DataFrames
    len_list = []

    # Iterate over samples in the specified path
    for samp in os.listdir(amr_annotations.path):
        anno_txt = os.path.join(amr_annotations.path, samp, "allele_mapping_data.txt")

        # Read each DataFrame and append it to the list
        len_sample = pd.read_csv(
            anno_txt, sep="\t", usecols=["Reference Sequence", "Reference Length"]
        )
        # index_col="Reference Sequence")

        len_list.append(len_sample)

    # Concatenate all DataFrames into a single DataFrame and remove duplicate rows
    len_df = pd.concat(len_list)
    len_df.drop_duplicates(inplace=True)
    len_df.columns = ["gene_id", "gene_length"]

    counts_df = table.to_dataframe()
    counts_df = counts_df.T

    with tempfile.TemporaryDirectory() as tmp:
        len = os.path.join(tmp, "len.csv")
        counts = os.path.join(tmp, "counts.csv")
        len_df.to_csv(len, index=False)
        counts_df.to_csv(counts, index=True)

        run_rnanorm_tpm(cwd=tmp, len=len, counts=counts)
        normalized_table = pd.read_csv(
            os.path.join(tmp, "out.csv"),
            index_col=0,
        )
        normalized_table.index.name = "sample_id"
    return normalized_table


def run_rnanorm_tpm(
    cwd: str,
    len: str,
    counts: str,
):
    cmd = ["rnanorm", "tpm", counts, "--gene-lengths", len, "--out", f"{cwd}/out.csv"]

    try:
        run_command(cmd, cwd, verbose=True, shell=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running RNAnorm, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
