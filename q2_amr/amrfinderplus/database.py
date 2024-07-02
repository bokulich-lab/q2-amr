import os
import subprocess

from qiime2.util import duplicate

from q2_amr.amrfinderplus.types import AMRFinderPlusDatabaseDirectoryFormat
from q2_amr.card.utils import run_command


def fetch_amrfinderplus_db() -> AMRFinderPlusDatabaseDirectoryFormat:
    amrfinderplus_db = AMRFinderPlusDatabaseDirectoryFormat()

    # Run AMRFinderPlus u function that downloads the database
    run_amrfinder_u()

    # Define path where the database will be downloaded to
    conda_prefix = os.getenv("CONDA_PREFIX")
    amrfinder_db_path = os.path.join(
        conda_prefix, "share", "amrfinderplus", "data", "latest"
    )

    # Copy all files from amrfinder_db_path to database directory format
    _copy_all(amrfinder_db_path, amrfinderplus_db.path)

    return amrfinderplus_db


def _copy_all(src_dir, des_dir):
    # Copies all files from source directory to destination directory
    for file in os.listdir(src_dir):
        src = os.path.join(src_dir, file)
        des = os.path.join(des_dir, file)
        duplicate(src, des)


def run_amrfinder_u():
    cmd = ["amrfinder", "-u"]

    try:
        run_command(cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running AMRFinderPlus, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
