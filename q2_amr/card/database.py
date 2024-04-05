import glob
import gzip
import os
import shutil
import subprocess
import tarfile
import tempfile

import requests
from tqdm import tqdm

from q2_amr.card.utils import run_command
from q2_amr.types._format import (
    CARDDatabaseDirectoryFormat,
    CARDKmerDatabaseDirectoryFormat,
)


def fetch_card_db() -> (CARDDatabaseDirectoryFormat, CARDKmerDatabaseDirectoryFormat):
    # Fetch CARD and WildCARD data from CARD website
    try:
        response_card = requests.get(
            "https://card.mcmaster.ca/latest/data", stream=True
        )
        response_wildcard = requests.get(
            "https://card.mcmaster.ca/latest/variants", stream=True
        )

        # Get content length to calculate progress
        total_size_card = int(response_card.headers.get("content-length", 0))
        total_size_wildcard = int(response_wildcard.headers.get("content-length", 0))

        # Initialize progress bars
        progress_bar_card = tqdm(
            total=total_size_card,
            unit="B",
            unit_scale=True,
            desc="Downloading CARD database",
        )

        progress_bar_wildcard = tqdm(
            total=total_size_wildcard,
            unit="B",
            unit_scale=True,
            desc="Downloading WildCARD database",
        )

        # Update progress bars

        for data_card in response_card.iter_content(chunk_size=4096):
            progress_bar_card.update(len(data_card))
        progress_bar_card.close()

        for data_wildcard in response_wildcard.iter_content(chunk_size=4096):
            progress_bar_wildcard.update(len(data_wildcard))
        progress_bar_wildcard.close()

    except requests.ConnectionError as e:
        error_msg = "Unable to connect to the CARD server. Please try again later."
        raise requests.ConnectionError(error_msg) from e

    msg = "Download finished. Extracting database files..."
    print(f"{msg}", end="", flush=True)

    # Create temporary directory for WildCARD data
    with tempfile.TemporaryDirectory() as tmp_dir:
        os.mkdir(os.path.join(tmp_dir, "wildcard"))

        # Extract tar.bz2 archives and store files in dirs "card" and "wildcard_zip"
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

        # List of files to be unzipped
        files = (
            "index-for-model-sequences.txt.gz",
            "nucleotide_fasta_protein_homolog_model_variants.fasta.gz",
            "nucleotide_fasta_protein_overexpression_model_variants.fasta.gz",
            "nucleotide_fasta_protein_variant_model_variants.fasta.gz",
            "nucleotide_fasta_rRNA_gene_variant_model_variants.fasta.gz",
            "61_kmer_db.json.gz",
            "all_amr_61mers.txt.gz",
        )

        # Unzip gzip files and save them in "wildcard" dir
        for file in files:
            with gzip.open(
                os.path.join(tmp_dir, "wildcard_zip", file), "rb"
            ) as f_in, open(
                os.path.join(tmp_dir, "wildcard", file[:-3]), "wb"
            ) as f_out:
                f_out.write(f_in.read())

        print(f"\r{msg} Done.", flush=True)

        msg = "Preprocessing database files..."
        print(f"{msg}", end="", flush=True)

        # Preprocess data for CARD and WildCARD
        # This creates additional fasta files in the temp directory
        preprocess(dir=tmp_dir, operation="card")
        preprocess(dir=tmp_dir, operation="wildcard")

        print(f"\r{msg} Done.", flush=True)

        msg = "Creating database artifacts..."
        print(f"{msg}", end="", flush=True)

        # Create CARD and Kmer database objects
        card_db = CARDDatabaseDirectoryFormat()
        kmer_db = CARDKmerDatabaseDirectoryFormat()

        # Find names of CARD database files created by preprocess function
        card_db_files = [
            os.path.basename(file)
            for file in glob.glob(os.path.join(tmp_dir, "card_database_v*.fasta"))
        ]

        # Lists of filenames to be moved to CARD and Kmer database objects
        wildcard_to_card_db = [
            "index-for-model-sequences.txt",
            "nucleotide_fasta_protein_homolog_model_variants.fasta",
            "nucleotide_fasta_protein_overexpression_model_variants.fasta",
            "nucleotide_fasta_protein_variant_model_variants.fasta",
            "nucleotide_fasta_rRNA_gene_variant_model_variants.fasta",
        ]
        tmp_to_card_db = [
            "wildcard_database_v0.fasta",
            "wildcard_database_v0_all.fasta",
            card_db_files[0],
            card_db_files[1],
        ]
        wildcard_to_kmer_db = ["all_amr_61mers.txt", "61_kmer_db.json"]

        # List of source and destination paths for files
        src_des_list = [
            (os.path.join(tmp_dir, "card"), str(card_db)),
            (os.path.join(tmp_dir, "wildcard"), str(card_db)),
            (tmp_dir, str(card_db)),
            (os.path.join(tmp_dir, "wildcard"), str(kmer_db)),
        ]

        # Move all files from source path to destination path
        for file_list, src_des in zip(
            [["card.json"], wildcard_to_card_db, tmp_to_card_db, wildcard_to_kmer_db],
            src_des_list,
        ):
            for file in file_list:
                shutil.move(
                    os.path.join(src_des[0], file), os.path.join(src_des[1], file)
                )

        return card_db, kmer_db


def preprocess(dir, operation):
    if operation == "card":
        # Run RGI command for CARD data
        cmd = ["rgi", "card_annotation", "-i", "card/card.json"]
    elif operation == "wildcard":
        # Run RGI command for WildCARD data
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
