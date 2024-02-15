import glob
import os
import shutil
import subprocess
import tempfile

from q2_amr.card.utils import load_card_db, run_command
from q2_amr.types import (
    CARDAnnotationDirectoryFormat,
    CARDDatabaseDirectoryFormat,
    CARDKmerDatabaseDirectoryFormat,
)
from q2_amr.types._format import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
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
    amr_annotations: CARDAlleleAnnotationDirectoryFormat,
    card_db: CARDDatabaseDirectoryFormat,
    kmer_db: CARDKmerDatabaseDirectoryFormat,
    minimum: int = 10,
    threads: int = 1,
) -> (
    CARDReadsAlleleKmerAnalysisDirectoryFormat,
    CARDReadsGeneKmerAnalysisDirectoryFormat,
):
    kmer_analysis = kmer_query(
        "reads", card_db, kmer_db, amr_annotations, minimum, threads
    )
    return kmer_analysis


def kmer_query(data_type, card_db, kmer_db, amr_annotations, minimum, threads):
    if data_type == "reads":
        annotation_file = "sorted.length_100.bam"
        input_type = "bwt"
        reads_allele_kmer_analysis = CARDReadsAlleleKmerAnalysisDirectoryFormat()
        reads_gene_kmer_analysis = CARDReadsGeneKmerAnalysisDirectoryFormat()
        kmer_analysis = (reads_allele_kmer_analysis, reads_gene_kmer_analysis)
    else:
        annotation_file = "amr_annotation.json"
        input_type = "rgi"
        mags_kmer_analysis = CARDMAGsKmerAnalysisDirectoryFormat()
        kmer_analysis = mags_kmer_analysis
    with tempfile.TemporaryDirectory() as tmp:
        # Load all necessary files and retrieve K-mer size
        load_card_db(tmp=tmp, card_db=card_db, kmer_db=kmer_db, kmer=True)
        path_kmer_json = glob.glob(os.path.join(str(kmer_db), "*_kmer_db.json"))[0]
        kmer_size = os.path.basename(path_kmer_json).split("_")[0]

        # Retrieve all paths to annotation files
        for root, dirs, files in os.walk(str(amr_annotations)):
            if annotation_file in files:
                input_path = os.path.join(root, annotation_file)
                # Run kmer_query
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
                    files = (
                        f"output_{kmer_size}mer_analysis.allele.txt",
                        f"output_{kmer_size}mer_analysis.json",
                        f"output_{kmer_size}mer_analysis.gene.txt",
                    )
                    des_dir_allele = os.path.join(
                        str(reads_allele_kmer_analysis), path_split[-2]
                    )
                    des_dir_gene = os.path.join(
                        str(reads_gene_kmer_analysis), path_split[-2]
                    )
                    des_dirs = [des_dir_allele, des_dir_allele, des_dir_gene]

                else:
                    files = (
                        f"output_{kmer_size}mer_analysis_{input_type}_summary.txt",
                        f"output_{kmer_size}mer_analysis.json",
                    )
                    des_dir = os.path.join(
                        str(mags_kmer_analysis), path_split[-3], path_split[-2]
                    )
                    des_dirs = [des_dir, des_dir]

                for file, des_dir in zip(files, des_dirs):
                    os.makedirs(des_dir, exist_ok=True)
                    des_path = os.path.join(des_dir, file[7:])
                    shutil.move(os.path.join(tmp, file), des_path)
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
        str(minimum),
        "--threads",
        str(threads),
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
