import glob
import json
import os
import subprocess

from qiime2.util import duplicate

from q2_amr.utils.utils import run_command


def load_card_db(
    card_db,
    kmer_db=None,
    kmer: bool = False,
    fasta: bool = False,
    include_other_models: bool = False,
    include_wildcard: bool = False,
):
    # Get path to card.json
    path_card_json = str(card_db.path / "card.json")

    # Base command that only loads card.json
    cmd = ["rgi", "load", "--card_json", path_card_json]

    # Define suffixes for card fasta file
    models = ("_all", "_all_models") if include_other_models is True else ("", "")

    # Extend base command with flag to load card fasta file
    if fasta:
        # Retrieve the database version number from card.jason file
        with open(path_card_json) as f:
            card_data = json.load(f)
            version = card_data["_version"]

        # Define path to card fasta file
        path_card_fasta = os.path.join(
            str(card_db), f"card_database_v{version}{models[0]}.fasta"
        )

        # Extend base command
        cmd.extend([f"--card_annotation{models[1]}", path_card_fasta])

    # Extend base command with flag to load wildcard fasta file and index
    if include_wildcard:
        cmd.extend(
            [
                f"--wildcard_annotation{models[1]}",
                os.path.join(str(card_db), f"wildcard_database_v0{models[0]}.fasta"),
                "--wildcard_index",
                os.path.join(str(card_db), "index-for-model-sequences.txt"),
            ]
        )
    # Extend base command with flag to load kmer json and txt database files
    kmer_size = None
    if kmer:
        path_kmer_json = glob.glob(os.path.join(str(kmer_db), "*_kmer_db.json"))[0]
        path_kmer_txt = glob.glob(os.path.join(str(kmer_db), "all_amr_*mers.txt"))[0]
        kmer_size = os.path.basename(path_kmer_json).split("_")[0]
        cmd.extend(
            [
                "--kmer_database",
                path_kmer_json,
                "--amr_kmers",
                path_kmer_txt,
                "--kmer_size",
                os.path.basename(path_kmer_json).split("_")[0],
            ]
        )

    # Run command
    try:
        run_command(cmd=cmd, cwd=None, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
    return kmer_size


def copy_files(file_paths: list, *dst_path_components):
    """
    Creates a destination file path out of the *dst_path_components. Then creates
    the directory for the destination file path if it doesn't exist already and
    finally copies the file from source path to destination path.

    Args:
        file_paths (list): A list of source file paths to be copied.
        *dst_path_components: Variable number of arguments representing destination
        path components that will be joined together to form the destination file
        path.
    """
    for src in file_paths:
        # Construct destination file path with destination file path components
        dst = os.path.join(*dst_path_components, os.path.basename(src))

        # Create destination directory if it not already exists
        os.makedirs(os.path.dirname(dst), exist_ok=True)

        # Copy file from source to destination
        duplicate(src, dst)
