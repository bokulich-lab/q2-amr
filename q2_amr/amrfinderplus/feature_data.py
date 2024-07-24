import glob
import os
import shutil
import tempfile

from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.genome_data import (
    GenesDirectoryFormat,
    LociDirectoryFormat,
    ProteinsDirectoryFormat,
)

from q2_amr.amrfinderplus.types import (
    AMRFinderPlusAnnotationsDirFmt,
    AMRFinderPlusDatabaseDirFmt,
)
from q2_amr.amrfinderplus.utils import run_amrfinderplus_n


def _validate_inputs(mags, loci, proteins):
    if mags and loci and not proteins:
        raise ValueError(
            "Loci input can only be given in combination with proteins input."
        )
    if mags and not loci and proteins:
        raise ValueError(
            "MAGs and proteins inputs together can only "
            "be given in combination with loci input."
        )
    if not mags and not proteins:
        raise ValueError("MAGs or proteins input has to be provided.")


def _get_file_paths(file, mags, proteins, loci):
    if mags:
        mag_id = os.path.splitext(os.path.basename(file))[0]
        if proteins:
            protein_path = os.path.join(str(proteins), f"{mag_id}_proteins.fasta")
            if not os.path.exists(protein_path):
                raise ValueError(
                    f"Proteins file for ID '{mag_id}' is missing in proteins input."
                )
        else:
            protein_path = None
    elif proteins:
        mag_id = os.path.splitext(os.path.basename(file))[0][:-9]
        protein_path = file
    else:
        raise ValueError("Either mags or proteins must be provided.")

    if loci:
        gff_path = os.path.join(str(loci), f"{mag_id}_loci.gff")
        if not os.path.exists(gff_path):
            raise ValueError(f"GFF file for ID '{mag_id}' is missing in loci input.")
    else:
        gff_path = None

    return mag_id, protein_path, gff_path


def _move_or_create_files(src_dir: str, mag_id: str, file_operations: dict):
    for file_name, target_dir in file_operations:
        if os.path.exists(os.path.join(src_dir, file_name)):
            shutil.move(
                os.path.join(src_dir, file_name),
                os.path.join(str(target_dir), f"{mag_id}_{file_name}"),
            )
        else:
            with open(os.path.join(str(target_dir), f"{mag_id}_{file_name}"), "w"):
                pass


def annotate_feature_data_amrfinderplus(
    amrfinderplus_db: AMRFinderPlusDatabaseDirFmt,
    mags: MAGSequencesDirFmt = None,
    proteins: ProteinsDirectoryFormat = None,
    loci: LociDirectoryFormat = None,
    organism: str = None,
    plus: bool = False,
    report_all_equal: bool = False,
    ident_min: float = None,
    curated_ident: bool = False,
    coverage_min: float = 0.5,
    translation_table: str = "11",
    annotation_format: str = "prodigal",
    report_common: bool = False,
    gpipe_org: bool = False,
    threads: int = None,
) -> (
    AMRFinderPlusAnnotationsDirFmt,
    AMRFinderPlusAnnotationsDirFmt,
    GenesDirectoryFormat,
    ProteinsDirectoryFormat,
):
    # Check for unallowed input combinations
    _validate_inputs(mags, loci, proteins)

    # Create all output directory formats
    amr_annotations = AMRFinderPlusAnnotationsDirFmt()
    amr_all_mutations = AMRFinderPlusAnnotationsDirFmt()
    amr_genes = GenesDirectoryFormat()
    amr_proteins = ProteinsDirectoryFormat()

    if mags:
        files = glob.glob(os.path.join(str(mags), "*"))
    elif proteins:
        files = glob.glob(os.path.join(str(proteins), "*"))

    with tempfile.TemporaryDirectory() as tmp:
        for file in files:
            mag_id, protein_path, gff_path = _get_file_paths(file, mags, proteins, loci)

            # Run amrfinderplus
            run_amrfinderplus_n(
                working_dir=tmp,
                amrfinderplus_db=amrfinderplus_db,
                dna_sequences=file if mags else None,
                protein_sequences=protein_path,
                gff=gff_path,
                organism=organism,
                plus=plus,
                report_all_equal=report_all_equal,
                ident_min=ident_min,
                curated_ident=curated_ident,
                coverage_min=coverage_min,
                translation_table=translation_table,
                annotation_format=annotation_format,
                report_common=report_common,
                gpipe_org=gpipe_org,
                threads=threads,
            )

            # Move mutations, genes and proteins files from tmp dir to the output
            # directory format, if organism, dna_sequence and proteins parameters
            # are specified. Else create empty placeholder files.
            file_operations = [
                ("amr_annotations.tsv", amr_annotations),
                ("amr_all_mutations.tsv", amr_all_mutations),
                ("amr_genes.fasta", amr_genes),
                ("amr_proteins.fasta", amr_proteins),
            ]

            # Loop through each file operation
            _move_or_create_files(tmp, mag_id, file_operations)

    return amr_annotations, amr_all_mutations, amr_genes, amr_proteins
