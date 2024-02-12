import os

import biom
import pandas as pd
from rnanorm import CPM, CTF, CUF, FPKM, TMM, TPM, UQ

from q2_amr.card.utils import InvalidParameterCombinationError
from q2_amr.types import GeneLengthDirectoryFormat


def normalize(
    table: biom.Table,
    method: str,
    m_trim: float = 0.3,
    a_trim: float = 0.05,
    gene_length: GeneLengthDirectoryFormat = None,
) -> pd.DataFrame:
    # Create Dataframe with counts from biom.Table
    counts = pd.DataFrame(
        data=table.matrix_data.toarray(),
        index=table.ids(axis="observation"),
        columns=table.ids(axis="sample"),
    ).T
    if method in ["tpm", "fpkm", "uq", "cuf", "cpm"]:
        # Raise Error if m or a-trim parameters are given with methods TPM, FPKM, UQ,
        # CPM or CUF
        if m_trim != 0.3 or a_trim != 0.05:
            raise InvalidParameterCombinationError(
                "Parameters m-trim and a-trim can only be used with methods TMM and "
                "CTF."
            )
        if method in ["tpm", "fpkm"]:
            # Raise Error if gene-length is missing when using methods TPM or FPKM
            if not gene_length:
                raise ValueError("gene-length input is missing.")
            # Create pd.Series from gene_length input
            gene_length_series = pd.read_csv(
                os.path.join(gene_length.path, "gene_length.txt"),
                sep="\t",
                header=None,
                names=["index", "values"],
                index_col="index",
                squeeze=True,
            )
            # Raise Error if there are genes in the counts that are not present in the
            # gene length
            if not set(counts.columns).issubset(set(gene_length_series.index)):
                only_in_counts = set(counts.columns) - set(gene_length_series.index)
                raise ValueError(
                    f"There are genes present in the FeatureTable that are not present "
                    f"in the gene-length input. Missing lengths for genes: "
                    f"{only_in_counts}"
                )
            # Define the methods TPM and FPKM with the gene length series as an input
            methods = {
                "tpm": TPM(gene_lengths=gene_length_series),
                "fpkm": FPKM(gene_lengths=gene_length_series),
            }
    if method in ["tmm", "uq", "cuf", "ctf", "cpm"]:
        # Raise Error if gene-length is given when using methods TMM, UQ, CUF, CPM or
        # CTF
        if gene_length:
            raise ValueError(
                "gene-length input can only be used with FPKM and TPM methods."
            )
        # Define the methods TMM and CTF with parameters, also UQ, CPM and CUF
        methods = {
            "tmm": TMM(m_trim=m_trim, a_trim=a_trim),
            "ctf": CTF(m_trim=m_trim, a_trim=a_trim),
            "uq": UQ(),
            "cuf": CUF(),
            "cpm": CPM(),
        }
    # Run normalization method on count dataframe
    normalized = methods[method].set_output(transform="pandas").fit_transform(counts)
    normalized.index.name = "sample_id"
    return normalized
