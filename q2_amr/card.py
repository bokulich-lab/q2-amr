import json
import os
import shutil
import subprocess
import tarfile
import tempfile
from distutils.dir_util import copy_tree

import pandas as pd
import pkg_resources
import q2templates
import requests
import skbio
from q2_types.feature_data import ProteinFASTAFormat, DNAFASTAFormat
from qiime2.core.type import List

from skbio import Protein, DNA

from q2_amr.types import CARDAnnotationJSONFormat, CARDDatabaseFormat, CARDAnnotationTXTFormat, \
    CARDAnnotationDirectoryFormat
from q2_amr.utils import run_command

CARD_URL = "https://card.mcmaster.ca/download/0/broadstreet-v{}.tar.bz2"


def fetch_card_db(version: str = '3.2.6') -> pd.DataFrame:
    url = CARD_URL.format(version)
    try:
        response = requests.get(url, stream=True)
    except requests.ConnectionError as e:
        raise requests.ConnectionError('Network connectivity problems.') from e
    with tempfile.TemporaryDirectory() as tmp_dir:
        try:
            with tarfile.open(fileobj=response.raw, mode="r|bz2") as tar:
                tar.extractall(path=tmp_dir)
        except tarfile.ReadError as a:
            raise tarfile.ReadError('Tarfile is invalid.') from a
        card_path = os.path.join(tmp_dir, 'card.json')
        card_df = pd.read_json(card_path).transpose()
        return card_df


def annotate_card(input_sequence: DNAFASTAFormat,
                  alignment_tool: str = 'BLAST',
                  input_type: str = 'contig',
                  split_prodigal_jobs: bool = False,
                  include_loose: bool = False,
                  exclude_nudge: bool = False,
                  low_quality: bool = False,
                  num_threads: int = 8) -> (CARDAnnotationDirectoryFormat, ProteinFASTAFormat, DNAFASTAFormat):
    with tempfile.TemporaryDirectory() as tmp:
        run_rgi_main(tmp, input_sequence, alignment_tool, input_type, split_prodigal_jobs, include_loose, exclude_nudge,
                     low_quality, num_threads)
        amr_annotation_df = pd.read_csv(f'{tmp}/output.txt', sep="\t")
        amr_annotations = CARDAnnotationDirectoryFormat()
        shutil.move(f'{tmp}/output.txt', f"{str(amr_annotations)}/amr_annotation.txt")
        shutil.move(f'{tmp}/output.json', f"{str(amr_annotations)}/amr_annotation.json")
    protein_annotation, dna_annotation = card_annotation_df_to_fasta(amr_annotation_df)
    return amr_annotations, protein_annotation, dna_annotation


def run_rgi_main(tmp,
                 input_sequence: DNAFASTAFormat,
                 alignment_tool: str = 'BLAST',
                 input_type: str = 'contig',
                 split_prodigal_jobs: bool = False,
                 include_loose: bool = False,
                 exclude_nudge: bool = False,
                 low_quality: bool = False,
                 num_threads: int = 8):
    cmd = ['rgi', 'main', '--input_sequence', f'{str(input_sequence)}', '--output_file', f'{tmp}/output', '-n',
           f'{num_threads}']
    if include_loose:
        cmd.append("--include_loose")
    if not exclude_nudge:
        cmd.append("--exclude_nudge")
    if low_quality:
        cmd.append("--low_quality")
    if split_prodigal_jobs:
        cmd.append("--split_prodigal_jobs")
    cmd.extend(['--alignment_tool', f'{alignment_tool}'])
    cmd.extend(['--input_type', f'{input_type}'])
    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def card_annotation_df_to_fasta(input_df: pd.DataFrame):
    protein_fasta = ProteinFASTAFormat()
    dna_fasta = DNAFASTAFormat()
    with open(str(protein_fasta), 'a') as proteinf, open(str(dna_fasta), 'a') as dnaf:
        for index, row in input_df.iterrows():
            protein_object = Protein(row['Predicted_Protein'])
            protein_object.metadata['id'] = row['ORF_ID']
            protein_object.metadata['description'] = row['ARO']
            skbio.io.write(protein_object, format='fasta', into=proteinf)
            dna_object = DNA(row['Predicted_DNA'])
            dna_object.metadata['id'] = row['ORF_ID']
            dna_object.metadata['description'] = row['ARO']
            skbio.io.write(dna_object, format='fasta', into=dnaf)
    return protein_fasta, dna_fasta


def heatmap(output_dir: str,
            amr_annotation: CARDAnnotationJSONFormat,
            clus: str = None,
            cat: str = None,
            display: str = 'plain',
            frequency: bool = False):
    TEMPLATES = pkg_resources.resource_filename("q2_amr", "assets")
    with tempfile.TemporaryDirectory() as tmp:
        results_dir = os.path.join(tmp, "results")
        json_files_dict = os.path.join(tmp, "json_files")
        os.makedirs(results_dir)
        os.makedirs(json_files_dict)
        for file in amr_annotation:
            shutil.copy(file, json_files_dict)
        run_rgi_heatmap(tmp, json_files_dict, clus, cat, display, frequency)
        change_names(results_dir)
        copy_tree(os.path.join(TEMPLATES, "rgi"), output_dir)
        copy_tree(results_dir, os.path.join(output_dir, "rgi_data"))
    context = {
        "tabs": [
            {"title": "QC report", "url": "index.html"}]
    }
    index = os.path.join(TEMPLATES, 'rgi', 'index.html')
    templates = [index]
    q2templates.render(templates, output_dir, context=context)


def run_rgi_heatmap(tmp, json_files_dict, clus, cat, display, frequency):
    cmd = ['rgi', 'heatmap', '--input', json_files_dict, '--output', f'{tmp}/results/heatmap',
           "--display", str(display)]
    if clus:
        cmd.extend(["--clus", str(clus)])
    if cat:
        cmd.extend(["--cat", str(cat)])
    if frequency:
        cmd.append("--frequency")
    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def change_names(results_dir):
    extensions = [".eps", ".csv", ".png"]  # Replace with the file extensions you want to rename
    # Get all the files in the directory
    files = os.listdir(results_dir)
    # Loop through all the files in the directory
    for filename in files:
        # Check if the file extension is in the list of extensions you want to rename
        if os.path.splitext(filename)[1] in extensions:
            # Construct the new file name with the correct extension
            file_ext = os.path.splitext(filename)[1]
            new_filename = "heatmap" + file_ext
            # Rename the file
            old_path = os.path.join(results_dir, filename)
            new_path = os.path.join(results_dir, new_filename)
            os.rename(old_path, new_path)
