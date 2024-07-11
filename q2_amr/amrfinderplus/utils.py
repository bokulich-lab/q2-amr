import subprocess

from q2_amr.card.utils import run_command


def run_amrfinderplus_n(
    working_dir,
    amrfinderplus_db,
    dna_sequence,
    protein_sequence,
    gff,
    organism,
    plus,
    report_all_equal,
    ident_min,
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
    if dna_sequence:
        cmd.extend(
            [
                "-n",
                dna_sequence,
                "--nucleotide_output",
                f"{working_dir}/amr_genes.fasta",
            ]
        )
    if protein_sequence:
        cmd.extend(
            [
                "-p",
                protein_sequence,
                "--protein_output",
                f"{working_dir}/amr_proteins.fasta",
            ]
        )
    if gff:
        cmd.extend(["-g", gff])
    if threads:
        cmd.extend(["--threads", str(threads)])
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
    if ident_min:
        cmd.extend(["--ident_min", str(ident_min)])
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
