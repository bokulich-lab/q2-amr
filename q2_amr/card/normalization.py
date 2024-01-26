import os

import biom
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from rnanorm import TPM

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
    # Initialize an empty series for gene lengths
    len_all = pd.Series()

    # Iterate over samples in the specified path
    for samp in os.listdir(amr_annotations.path):
        anno_txt = os.path.join(amr_annotations.path, samp, "allele_mapping_data.txt")

        # Read each DataFrame and append it to the list
        len_sample = pd.read_csv(
            anno_txt, sep="\t", usecols=["Reference Sequence", "Reference Length"]
        ).set_index("Reference Sequence")["Reference Length"]

        len_all = len_all.combine_first(len_sample)
    len_all = len_all
    # len_all = len_all.values.reshape(-1, 1)
    counts_arr = table.matrix_data.toarray().T
    # counts_arr = counts_arr.T
    # len_all = len_all.reshape(-1, 1)
    tpm_normalizer = TPM(gene_lengths=len_all).set_output(transform="pandas")
    tpm_df = tpm_normalizer.fit_transform(counts_arr)
    tpm_df.index.name = "sample_id"
    return tpm_df
