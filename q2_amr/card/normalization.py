import biom
import pandas as pd
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
    # dds = DeseqDataSet(counts=df, metadata=metadata)
    # dds.fit_size_factors()

    # normalizedcounts = dds.obsm["size_factors"]
    return df


def normalize_tpm(
    table: biom.Table, amr_annotations: CARDAlleleAnnotationDirectoryFormat
) -> pd.DataFrame:
    # # Initialize an empty series for gene lengths
    # len_all = pd.Series()
    #
    # # Iterate over samples in the specified path
    # for samp in os.listdir(amr_annotations.path):
    #     anno_txt = os.path.join(amr_annotations.path, samp, "allele_mapping_data.txt")
    #
    #     # Read each DataFrame and append it to the list
    #     len_sample = pd.read_csv(
    #         anno_txt, sep="\t", usecols=["Reference Sequence", "Reference Length"]
    #     ).set_index("Reference Sequence")["Reference Length"]
    #
    #     len_all = len_all.combine_first(len_sample)
    len_all = pd.read_csv(
        "/Users/rischv/Desktop/genelengts.txt",
        sep="\t",
        header=None,
        names=["index", "values"],
        index_col="index",
        squeeze=True,
    )
    len_all = len_all.astype(int)
    df = pd.DataFrame(
        data=table.matrix_data.toarray(),
        index=table.ids(axis="observation"),
        columns=table.ids(axis="sample"),
    ).T
    # Work with pandas.
    # len_all.to_csv("/Users/rischv/Desktop/genelengts.txt", sep='\t', index=True,
    # header=False)

    transformer = TPM(gene_lengths=len_all).set_output(transform="pandas")
    exp_normalized = transformer.fit_transform(df)

    exp_normalized.index.name = "sample_id"
    return exp_normalized
