import os

import pandas as pd
from q2_types.feature_data import SequenceCharacteristicsDirectoryFormat
from rnanorm import CPM, CTF, CUF, FPKM, TMM, TPM, UQ


def normalize(
    table: pd.DataFrame,
    method: str,
    m_trim: float = None,
    a_trim: float = None,
    gene_length: SequenceCharacteristicsDirectoryFormat = None,
) -> pd.DataFrame:
    # Validate parameter combinations and set trim parameters
    m_trim, a_trim = _validate_parameters(method, m_trim, a_trim, gene_length)

    # Process gene_lengths input and define methods that need gene_lengths input
    if method in ["tpm", "fpkm"]:
        lengths = _convert_lengths(table, gene_length)

        methods = {
            "tpm": TPM(gene_lengths=lengths),
            "fpkm": FPKM(gene_lengths=lengths),
        }

    # Define remaining methods that don't need gene_lengths input
    else:
        methods = {
            "tmm": TMM(m_trim=m_trim, a_trim=a_trim),
            "ctf": CTF(m_trim=m_trim, a_trim=a_trim),
            "uq": UQ(),
            "cuf": CUF(),
            "cpm": CPM(),
        }

    # Run normalization method on frequency table
    normalized = methods[method].set_output(transform="pandas").fit_transform(table)
    normalized.index.name = "sample_id"

    return normalized


def _validate_parameters(method, m_trim, a_trim, gene_length):
    # Raise Error if gene-length is missing when using methods TPM or FPKM
    if method in ["tpm", "fpkm"] and not gene_length:
        raise ValueError("gene-length input is missing.")

    # Raise Error if gene-length is given when using methods TMM, UQ, CUF, CPM or CTF
    if method in ["tmm", "uq", "cuf", "ctf", "cpm"] and gene_length:
        raise ValueError(
            "gene-length input can only be used with FPKM and TPM methods."
        )

    # Raise Error if m_trim or a_trim are given when not using methods TMM or CTF
    if (method not in ["tmm", "ctf"]) and (m_trim is not None or a_trim is not None):
        raise ValueError(
            "Parameters m-trim and a-trim can only be used with methods TMM and CTF."
        )

    # Set m_trim and a_trim to their default values for methods TMM and CTF
    if method in ["tmm", "ctf"]:
        m_trim = 0.3 if m_trim is None else m_trim
        a_trim = 0.05 if a_trim is None else a_trim

    return m_trim, a_trim


def _convert_lengths(table, gene_length):
    # Read in table from sequence_characteristics.tsv as a pd.Series
    lengths = pd.read_csv(
        os.path.join(gene_length.path, "sequence_characteristics.tsv"),
        sep="\t",
        header=None,
        names=["index", "values"],
        index_col="index",
        squeeze=True,
        skiprows=1,
    )

    # Check if all gene IDs that are present in the table are also present in
    # the lengths
    if not set(table.columns).issubset(set(lengths.index)):
        only_in_counts = set(table.columns) - set(lengths.index)
        raise ValueError(
            f"There are genes present in the FeatureTable that are not present "
            f"in the gene-length input. Missing lengths for genes: "
            f"{only_in_counts}"
        )
    return lengths
