import glob
import gzip
import os
import shutil
import subprocess
import tarfile
import tempfile

import requests

from q2_amr.card.utils import run_command
from q2_amr.types._format import (
    CARDDatabaseDirectoryFormat,
    CARDKmerDatabaseDirectoryFormat,
)

CARD_URL = "https://card.mcmaster.ca/latest/data"
WILDCARD_URL = "https://card.mcmaster.ca/latest/variants"


def fetch_card_db() -> (CARDDatabaseDirectoryFormat, CARDKmerDatabaseDirectoryFormat):
    try:
        response_card = requests.get(CARD_URL, stream=True)
        response_wildcard = requests.get(WILDCARD_URL, stream=True)
    except requests.ConnectionError as e:
        raise requests.ConnectionError("Network connectivity problems.") from e
    with tempfile.TemporaryDirectory() as tmp_dir:
        os.mkdir(os.path.join(tmp_dir, "wildcard"))
        try:
            with tarfile.open(
                fileobj=response_card.raw, mode="r|bz2"
            ) as c_tar, tarfile.open(
                fileobj=response_wildcard.raw, mode="r|bz2"
            ) as wc_tar:
                c_tar.extractall(path=os.path.join(tmp_dir, "card"))
                wc_tar.extractall(path=os.path.join(tmp_dir, "wildcard_zip"))
        except tarfile.ReadError as a:
            raise tarfile.ReadError("Tarfile is invalid.") from a
        files = (
            "index-for-model-sequences.txt.gz",
            "nucleotide_fasta_protein_homolog_model_variants.fasta.gz",
            "nucleotide_fasta_protein_overexpression_model_variants.fasta.gz",
            "nucleotide_fasta_protein_variant_model_variants.fasta.gz",
            "nucleotide_fasta_rRNA_gene_variant_model_variants.fasta.gz",
            "61_kmer_db.json.gz",
            "all_amr_61mers.txt.gz",
        )
        for file in files:
            with gzip.open(
                os.path.join(tmp_dir, "wildcard_zip", file), "rb"
            ) as f_in, open(
                os.path.join(tmp_dir, "wildcard", file[:-3]), "wb"
            ) as f_out:
                f_out.write(f_in.read())

        preprocess(dir=tmp_dir, operation="card")
        preprocess(dir=tmp_dir, operation="wildcard")
        card_db = CARDDatabaseDirectoryFormat()
        kmer_db = CARDKmerDatabaseDirectoryFormat()
        card_db_files = [
            os.path.basename(file)
            for file in glob.glob(os.path.join(tmp_dir, "card_database_v*.fasta"))
        ]
        file_src_des = [
            ("card.json", os.path.join(tmp_dir, "card"), str(card_db)),
            (
                "index-for-model-sequences.txt",
                os.path.join(tmp_dir, "wildcard"),
                str(card_db),
            ),
            ("wildcard_database_v0.fasta", tmp_dir, str(card_db)),
            ("wildcard_database_v0_all.fasta", tmp_dir, str(card_db)),
            ("all_amr_61mers.txt", os.path.join(tmp_dir, "wildcard"), str(kmer_db)),
            ("61_kmer_db.json", os.path.join(tmp_dir, "wildcard"), str(kmer_db)),
            (
                "nucleotide_fasta_protein_homolog_model_variants.fasta",
                os.path.join(tmp_dir, "wildcard"),
                str(card_db),
            ),
            (
                "nucleotide_fasta_protein_overexpression_model_variants.fasta",
                os.path.join(tmp_dir, "wildcard"),
                str(card_db),
            ),
            (
                "nucleotide_fasta_protein_variant_model_variants.fasta",
                os.path.join(tmp_dir, "wildcard"),
                str(card_db),
            ),
            (
                "nucleotide_fasta_rRNA_gene_variant_model_variants.fasta",
                os.path.join(tmp_dir, "wildcard"),
                str(card_db),
            ),
            (card_db_files[0], tmp_dir, str(card_db)),
            (card_db_files[1], tmp_dir, str(card_db)),
        ]
        for file, src_dir, des_dir in file_src_des:
            shutil.move(os.path.join(src_dir, file), os.path.join(des_dir, file))
        return card_db, kmer_db


def preprocess(dir, operation):
    if operation == "card":
        cmd = ["rgi", "card_annotation", "-i", "card/card.json"]

    elif operation == "wildcard":
        cmd = [
            "rgi",
            "wildcard_annotation",
            "-i",
            "wildcard",
            "--card_json",
            "card/card.json",
            "-v",
            "0",
        ]
    try:
        run_command(cmd, dir, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
