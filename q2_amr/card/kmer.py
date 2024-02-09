import glob
import os
import shutil
import subprocess
import tempfile
from typing import Union

from q2_amr.card.utils import load_card_db, run_command
from q2_amr.types import (
    CARDAnnotationDirectoryFormat,
    CARDDatabaseDirectoryFormat,
    CARDKmerDatabaseDirectoryFormat,
)
from q2_amr.types._format import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDReadsKmerAnalysisDirectoryFormat,
)


def kmer_query_mags_card(
    amr_annotations: CARDAnnotationDirectoryFormat,
    kmer_db: CARDKmerDatabaseDirectoryFormat,
    card_db: CARDDatabaseDirectoryFormat,
    minimum: int = 10,
    threads: int = 1,
) -> CARDMAGsKmerAnalysisDirectoryFormat:
    kmer_analysis = kmer_query(
        "mags", card_db, kmer_db, amr_annotations, minimum, threads
    )
    return kmer_analysis


def kmer_query_reads_card(
    amr_annotations: Union[
        CARDAlleleAnnotationDirectoryFormat, CARDGeneAnnotationDirectoryFormat
    ],
    kmer_db: CARDKmerDatabaseDirectoryFormat,
    card_db: CARDDatabaseDirectoryFormat,
    minimum: int = 10,
    threads: int = 1,
) -> CARDReadsKmerAnalysisDirectoryFormat:
    kmer_analysis = kmer_query(
        "reads", card_db, kmer_db, amr_annotations, minimum, threads
    )
    return kmer_analysis


def kmer_query(data_type, card_db, kmer_db, amr_annotations, minimum, threads):
    if data_type == "reads":
        annotation_file = "sorted.length_100.bam"
        input_type = "bwt"
        kmer_analysis = CARDReadsKmerAnalysisDirectoryFormat()
    else:
        annotation_file = "amr_annotation.json"
        input_type = "rgi"
        kmer_analysis = CARDMAGsKmerAnalysisDirectoryFormat()
    with tempfile.TemporaryDirectory() as tmp:
        kmer_size = load_card_db(tmp=tmp, card_db=card_db, kmer_db=kmer_db, kmer=True)
        for root, dirs, files in os.walk(str(amr_annotations)):
            if annotation_file in files:
                input_path = os.path.join(root, annotation_file)
                run_rgi_kmer_query(
                    tmp=tmp,
                    input_file=input_path,
                    input_type=input_type,
                    kmer_size=kmer_size,
                    minimum=minimum,
                    threads=threads,
                )
                path_split = input_path.split(os.path.sep)
                if data_type == "reads":
                    des_dir = os.path.join(str(kmer_analysis), path_split[-2])
                else:
                    des_dir = os.path.join(
                        str(kmer_analysis), path_split[-3], path_split[-2]
                    )
                os.makedirs(des_dir)
                des_path = os.path.join(
                    des_dir, f"{kmer_size}mer_analysis_{data_type}.txt"
                )
                pattern = f"output_*mer_analysis_{input_type}_summary.txt"
                kmer_output = os.path.join(tmp, pattern)
                src_path = glob.glob(kmer_output)[0]
                shutil.move(src_path, des_path)
    return kmer_analysis


def run_rgi_kmer_query(tmp, input_file, input_type, kmer_size, minimum, threads):
    cmd = [
        "rgi",
        "kmer_query",
        "--input",
        input_file,
        f"--{input_type}",
        "--kmer_size",
        kmer_size,
        "--minimum",
        minimum,
        "--threads",
        threads,
        "--output",
        "output",
        "--local",
    ]

    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
