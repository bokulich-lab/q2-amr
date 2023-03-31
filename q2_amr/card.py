import os
import subprocess
import tarfile
import tempfile

import pandas as pd
import requests
import skbio
from q2_types.feature_data import ProteinFASTAFormat, DNAFASTAFormat
from skbio import Protein, DNA


def fetch_card_data(version: str = '3.2.6') -> pd.DataFrame:
    url = f"https://card.mcmaster.ca/download/0/broadstreet-v{version}.tar.bz2"
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


def card_annotation(sequences: DNAFASTAFormat,
                    alignment_tool: str = 'BLAST',
                    input_type: str = 'contig',
                    split_prodigal_jobs: bool = False,
                    loose: bool = False,
                    nudge: bool = True,
                    low_quality: bool = False,
                    threads: int = 8) -> (pd.DataFrame, ProteinFASTAFormat, DNAFASTAFormat):
    with tempfile.TemporaryDirectory() as tmp:
        cmd = [f'rgi main --input_sequence {str(sequences)} --output_file {tmp}/output -n {threads}']
        if loose:
            cmd.extend([" --include_loose"])
        if not nudge:
            cmd.extend([" --exclude_nudge"])
        if low_quality:
            cmd.extend([" --low_quality"])
        if split_prodigal_jobs:
            cmd.extend([" --split_prodigal_jobs"])
        if alignment_tool == 'Blast':
            cmd.extend([" -a BLAST"])
        elif alignment_tool == 'DIAMOND':
            cmd.extend([" -a DIAMOND"])
        if input_type == 'contig':
            cmd.extend([" -t contig"])
        elif input_type == 'protein':
            cmd.extend([" -t protein"])
        try:
            run_command(cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running rgi, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )
        rgi_output = pd.read_csv(f'{tmp}/output.txt', sep="\t")
    protein_fasta, dna_fasta = card_annotation_df_to_fasta(rgi_output)
    return rgi_output, protein_fasta, dna_fasta


def card_annotation_df_to_fasta(input_df):
    protein_fasta = ProteinFASTAFormat()
    dna_fasta = DNAFASTAFormat()
    with open(str(protein_fasta), 'a') as proteinf, open(str(dna_fasta), 'a') as dnaf:
        for index, row in input_df.iterrows():
            protein_object = Protein(row['Predicted_Protein'])
            dna_object = DNA(row['Predicted_DNA'])
            protein_object.metadata['id'] = row['ORF_ID']
            protein_object.metadata['description'] = row['ARO']
            dna_object.metadata['id'] = row['ORF_ID']
            dna_object.metadata['description'] = row['ARO']
            skbio.io.write(protein_object, format='fasta', into=proteinf)
            skbio.io.write(dna_object, format='fasta', into=dnaf)
    return protein_fasta, dna_fasta


def card_annotation_heatmap(rgi_input, output, clus: str = 'no', cat: str = 'no', frequency=False):
    cmd = [f'rgi heatmap --input {rgi_input} --output {output}']
    if frequency:
        cmd.extend(["--frequency"])
    if clus == 'both':
        cmd.extend(["-clus both"])
    elif clus == 'samples':
        cmd.extend(["-clus samples"])
    elif clus == 'genes':
        cmd.extend(["-clus genes"])
    elif clus == 'no':
        pass
    if cat == 'drug_class':
        cmd.extend(["-cat drug_class"])
    elif cat == 'resistance_mechanism':
        cmd.extend(["-cat resistance_mechanism"])
    elif cat == 'gene_family':
        cmd.extend(["-cat gene_family"])
    elif cat == 'no':
        pass
    try:
        run_command(cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


EXTERNAL_CMD_WARNING = (
    "Running external command line application(s). "
    "This may print messages to stdout and/or stderr.\n"
    "The command(s) being run are below. These commands "
    "cannot be manually re-run as they will depend on "
    "temporary files that no longer exist."
)


def run_command(cmd, verbose=True):
    if verbose:
        print(EXTERNAL_CMD_WARNING)
        print("\nCommand:", end=" ")
        print("".join(cmd), end="\n\n")
    subprocess.run(cmd, check=True, shell=True)
