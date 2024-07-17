import subprocess

from q2_amr.card.utils import run_command


def run_amrfinderplus_n(
    working_dir,
    amrfinderplus_db,
    dna_sequences,
    protein_sequences,
    gff,
    organism,
    plus,
    report_all_equal,
    ident_min,
    curated_ident,
    coverage_min,
    translation_table,
    threads,
):
    cmd = [
        "amrfinder",
        "--database",
        str(amrfinderplus_db),
        "-o",
        f"{working_dir}/amr_annotations.tsv",
        "--print_node",
    ]
    # Creates nucleotide fasta output if DNA sequences are given as input
    if dna_sequences:
        cmd.extend(
            [
                "-n",
                dna_sequences,
                "--nucleotide_output",
                f"{working_dir}/amr_genes.fasta",
            ]
        )
    # Creates protein fasta output if protein sequences are given as input
    if protein_sequences:
        cmd.extend(
            [
                "-p",
                protein_sequences,
                "--protein_output",
                f"{working_dir}/amr_proteins.fasta",
            ]
        )
    if gff:
        cmd.extend(["-g", gff])
    if threads:
        cmd.extend(["--threads", str(threads)])
    # Creates all mutations output if an organism is specified
    if organism:
        cmd.extend(
            [
                "--organism",
                organism,
                "--mutation_all",
                f"{working_dir}/amr_mutations.tsv",
            ]
        )
    if plus:
        cmd.append("--plus")
    if report_all_equal:
        cmd.append("--report_all_equal")
    # If curated_ident is True, it will overwrite the value specified with ident_min
    if ident_min and not curated_ident:
        cmd.extend(["--ident_min", str(ident_min)])
    if curated_ident:
        cmd.extend(["--ident_min", "-1"])
    if coverage_min:
        cmd.extend(["--coverage_min", str(coverage_min)])
    if translation_table:
        cmd.extend(["--translation_table", str(translation_table)])
    try:
        run_command(cmd, working_dir, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running AMRFinderPlus, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
