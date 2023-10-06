import glob
import gzip
import os
import shutil
import subprocess
import tarfile
import tempfile

import requests

from q2_amr.card.utils import run_command
from q2_amr.types import CARDDatabaseDirectoryFormat

CARD_URL = "https://card.mcmaster.ca/latest/data"
VARIANTS_URL = "https://card.mcmaster.ca/latest/variants"


def fetch_card_db() -> CARDDatabaseDirectoryFormat:
    try:
        response_card = requests.get(CARD_URL, stream=True)
        response_variants = requests.get(VARIANTS_URL, stream=True)
    except requests.ConnectionError as e:
        raise requests.ConnectionError("Network connectivity problems.") from e
    with tempfile.TemporaryDirectory() as tmp_dir:
        os.mkdir(os.path.join(tmp_dir, "variants_unzip"))
        try:
            with tarfile.open(fileobj=response_card.raw, mode="r|bz2") as tar:
                tar.extractall(path=os.path.join(tmp_dir, "card"))
        except tarfile.ReadError as a:
            raise tarfile.ReadError("Tarfile is invalid.") from a
        try:
            with tarfile.open(fileobj=response_variants.raw, mode="r|bz2") as tar:
                tar.extractall(path=os.path.join(tmp_dir, "variants"))
        except tarfile.ReadError as a:
            raise tarfile.ReadError("Tarfile is invalid.") from a

        # with tempfile.TemporaryDirectory() as tmp_dir:
        #     os.mkdir(os.path.join(tmp_dir, "variants_unzip"))
        #     shutil.copytree("/Users/vinzent/Desktop/bokulich_project/data/wildcard/card",
        #                     os.path.join(tmp_dir, "card"))
        #     shutil.copytree("/Users/vinzent/Desktop/bokulich_project/data/wildcard/variants",
        #                     os.path.join(tmp_dir, "variants"))
        files = (
            "index-for-model-sequences.txt.gz",
            "nucleotide_fasta_protein_homolog_model_variants.fasta.gz",
            "nucleotide_fasta_protein_overexpression_model_variants.fasta.gz",
            "nucleotide_fasta_protein_variant_model_variants.fasta.gz",
            "nucleotide_fasta_rRNA_gene_variant_model_variants.fasta.gz",
        )
        for file in files:
            with gzip.open(os.path.join(tmp_dir, "variants", file), "rb") as f_in, open(
                os.path.join(tmp_dir, "variants_unzip", file[:-3]), "wb"
            ) as f_out:
                f_out.write(f_in.read())
        preprocess(tmp_dir, "preprocess_card")
        preprocess(tmp_dir, "preprocess_variants")
        card_db = CARDDatabaseDirectoryFormat()
        shutil.move(
            os.path.join(tmp_dir, "card", "card.json"),
            os.path.join(str(card_db), "card.json"),
        )
        shutil.move(
            os.path.join(tmp_dir, "variants_unzip", "index-for-model-sequences.txt"),
            os.path.join(str(card_db), "index-for-model-sequences.txt"),
        )
        shutil.move(
            os.path.join(tmp_dir, "wildcard_database_v0.fasta"),
            os.path.join(str(card_db), "wildcard_database_v0.fasta"),
        )
        shutil.move(
            os.path.join(tmp_dir, "wildcard_database_v0_all.fasta"),
            os.path.join(str(card_db), "wildcard_database_v0_all.fasta"),
        )
        matching_files = glob.glob(os.path.join(tmp_dir, "card_database_v*.fasta"))
        for file_path in matching_files:
            shutil.move(
                file_path, os.path.join(str(card_db), os.path.basename(file_path))
            )
    return card_db


def preprocess(tmp, operation):
    if operation == "preprocess_card":
        cmd = ["rgi", "card_annotation", "-i", "card/card.json"]

    elif operation == "preprocess_variants":
        cmd = [
            "rgi",
            "wildcard_annotation",
            "-i",
            "variants_unzip",
            "--card_json",
            "card/card.json",
            "-v",
            "0",
        ]
    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
