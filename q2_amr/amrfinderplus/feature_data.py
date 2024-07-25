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
    # If mags is provided, mag_id is extracted from the file name.
    if mags:
        mag_id = os.path.splitext(os.path.basename(file))[0]

        # If proteins are provided, construct the expected protein file path.
        if proteins:
            protein_path = os.path.join(str(proteins), f"{mag_id}_proteins.fasta")

            # Raise an error if the expected protein file does not exist.
            if not os.path.exists(protein_path):
                raise ValueError(
                    f"Proteins file for ID '{mag_id}' is missing in proteins input."
                )
        else:
            protein_path = None

    # If only proteins are provided (without mags), determine mag_id and protein path.
    else:
        # Extract mag_id from the file name, excluding the last 9 characters
        # '_proteins'.
        mag_id = os.path.splitext(os.path.basename(file))[0][:-9]
        protein_path = file

    # If loci are provided, construct the expected GFF file path.
    if loci:
        gff_path = os.path.join(str(loci), f"{mag_id}_loci.gff")

        # Raise an error if the expected GFF file does not exist.
        if not os.path.exists(gff_path):
            raise ValueError(f"GFF file for ID '{mag_id}' is missing in loci input.")
    else:
        gff_path = None

    return mag_id, protein_path, gff_path


def _move_or_create_files(src_dir: str, mag_id: str, file_operations: list):
    # Loop through all files.
    for file_name, target_dir in file_operations:
        # If the file exists move it to the destination dir and attach mag_id.
        if os.path.exists(os.path.join(src_dir, file_name)):
            shutil.move(
                os.path.join(src_dir, file_name),
                os.path.join(str(target_dir), f"{mag_id}_{file_name}"),
            )
        # If the file does not exist, create empty placeholder file in the
        # destination dir.
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
    # Check for unallowed input combinations.
    _validate_inputs(mags, loci, proteins)

    # Create all output directories.
    amr_annotations = AMRFinderPlusAnnotationsDirFmt()
    amr_all_mutations = AMRFinderPlusAnnotationsDirFmt()
    amr_genes = GenesDirectoryFormat()
    amr_proteins = ProteinsDirectoryFormat()

    # Create list of files to loop over, if mags is provided then files in mags will be
    # used if only proteins is provided then files in proteins will be used
    if mags:
        files = glob.glob(os.path.join(str(mags), "*"))
    else:
        files = glob.glob(os.path.join(str(proteins), "*"))

    with tempfile.TemporaryDirectory() as tmp:
        # Loop over all files
        for file in files:
            # Get paths to protein and gff files, and get mag_id
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

            # Output filenames and output directories
            file_operations = [
                ("amr_annotations.tsv", amr_annotations),
                ("amr_all_mutations.tsv", amr_all_mutations),
                ("amr_genes.fasta", amr_genes),
                ("amr_proteins.fasta", amr_proteins),
            ]

            # Move the files or create placeholder files
            _move_or_create_files(tmp, mag_id, file_operations)

    return amr_annotations, amr_all_mutations, amr_genes, amr_proteins
