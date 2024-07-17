import os
import shutil
import tempfile

from q2_types.feature_data import (
    DNASequencesDirectoryFormat,
    ProteinSequencesDirectoryFormat,
)
from q2_types.genome_data import (
    GenesDirectoryFormat,
    LociDirectoryFormat,
    ProteinsDirectoryFormat,
)

from q2_amr.amrfinderplus.types import (
    AMRFinderPlusAnnotationDirFmt,
    AMRFinderPlusDatabaseDirFmt,
)
from q2_amr.amrfinderplus.utils import run_amrfinderplus_n


def annotate_sequences_amrfinderplus(
    amrfinderplus_db: AMRFinderPlusDatabaseDirFmt,
    dna_sequences: DNASequencesDirectoryFormat = None,
    protein_sequences: ProteinSequencesDirectoryFormat = None,
    gff: LociDirectoryFormat = None,
    organism: str = None,
    plus: bool = False,
    report_all_equal: bool = False,
    ident_min: float = None,
    curated_ident: bool = False,
    coverage_min: float = 0.5,
    translation_table: str = "11",
    threads: int = None,
) -> (
    AMRFinderPlusAnnotationDirFmt,
    AMRFinderPlusAnnotationDirFmt,
    GenesDirectoryFormat,
    ProteinsDirectoryFormat,
):
    # Check for unallowed input combinations
    if dna_sequences and gff and not protein_sequences:
        raise ValueError(
            "GFF input can only be given in combination with protein-sequence input."
        )
    if dna_sequences and not gff and protein_sequences:
        raise ValueError(
            "DNA-sequence and protein-sequence inputs together can only "
            "be given in combination with GFF input."
        )

    # Create all output directory formats
    annotations = AMRFinderPlusAnnotationDirFmt()
    mutations = AMRFinderPlusAnnotationDirFmt()
    genes = GenesDirectoryFormat()
    proteins = ProteinsDirectoryFormat()

    with tempfile.TemporaryDirectory() as tmp:
        # Run amrfinderplus function
        run_amrfinderplus_n(
            working_dir=tmp,
            amrfinderplus_db=amrfinderplus_db,
            dna_sequences=os.path.join(str(dna_sequences), "dna-sequences.fasta")
            if dna_sequences
            else None,
            protein_sequences=os.path.join(
                str(protein_sequences), "protein-sequences.fasta"
            )
            if protein_sequences
            else None,
            gff=os.path.join(str(gff), os.listdir(str(gff))[0]) if gff else None,
            organism=organism,
            plus=plus,
            report_all_equal=report_all_equal,
            ident_min=ident_min,
            curated_ident=curated_ident,
            coverage_min=coverage_min,
            translation_table=translation_table,
            threads=threads,
        )

        # Move annotations file from tmp dir to the output directory format
        shutil.move(os.path.join(tmp, "amr_annotations.tsv"), str(annotations))

        # Move mutations, genes and proteins files from tmp dir to the output
        # directory format, if organism, dna_sequence and protein_sequence parameters
        # are specified. Else create empty placeholder files.
        if organism:
            shutil.move(os.path.join(tmp, "amr_mutations.tsv"), str(mutations))
        else:
            with open(os.path.join(str(mutations), "amr_mutations.tsv"), "w"):
                pass

        if dna_sequences:
            shutil.move(os.path.join(tmp, "amr_genes.fasta"), str(genes))
        else:
            with open(os.path.join(str(genes), "amr_genes.fasta"), "w"):
                pass

        if protein_sequences:
            shutil.move(os.path.join(tmp, "amr_proteins.fasta"), str(proteins))
        else:
            with open(os.path.join(str(proteins), "amr_proteins.fasta"), "w"):
                pass

    return annotations, mutations, genes, proteins
