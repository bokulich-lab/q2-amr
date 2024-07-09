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
    AMRFinderPlusDatabaseDirFmt,
    ARMFinderPlusAnnotationDirFmt,
)
from q2_amr.amrfinderplus.utils import run_amrfinderplus_n


def annotate_sequences_amrfinderplus(
    amrfinderplus_db: AMRFinderPlusDatabaseDirFmt,
    dna_sequence: DNASequencesDirectoryFormat = None,
    protein_sequence: ProteinSequencesDirectoryFormat = None,
    gff: LociDirectoryFormat = None,
    organism: str = None,
    plus: bool = False,
    report_all_equal: bool = False,
    ident_min: float = None,
    coverage_min: float = 0.5,
    translation_table: str = "11",
    threads: int = None,
) -> (
    ARMFinderPlusAnnotationDirFmt,
    ARMFinderPlusAnnotationDirFmt,
    GenesDirectoryFormat,
    ProteinsDirectoryFormat,
):
    if dna_sequence and gff and not protein_sequence:
        raise ValueError(
            "GFF input can only be given in combination with " "protein-sequence input."
        )
    if dna_sequence and not gff and protein_sequence:
        raise ValueError(
            "DNA-sequence and protein-sequence inputs together can only "
            "be given in combination with GFF input."
        )

    annotations = ARMFinderPlusAnnotationDirFmt()
    mutations = ARMFinderPlusAnnotationDirFmt()
    genes = GenesDirectoryFormat()
    proteins = ProteinsDirectoryFormat()

    with tempfile.TemporaryDirectory() as tmp:
        run_amrfinderplus_n(
            working_dir=tmp,
            amrfinderplus_db=amrfinderplus_db,
            dna_sequence=os.path.join(str(dna_sequence), "dna-sequences.fasta")
            if dna_sequence
            else None,
            protein_sequence=os.path.join(
                str(protein_sequence), "protein-sequences.fasta"
            )
            if protein_sequence
            else None,
            gff=os.listdir(str(gff))[0] if gff else None,
            organism=organism,
            plus=plus,
            report_all_equal=report_all_equal,
            ident_min=ident_min,
            coverage_min=coverage_min,
            translation_table=translation_table,
            threads=threads,
            id="",
        )

        shutil.move(os.path.join(tmp, "amr_annotations.tsv"), str(annotations))

        if organism:
            shutil.move(os.path.join(tmp, "amr_mutations.tsv"), str(mutations))
        else:
            with open(os.path.join(str(mutations), "amr_mutations.tsv"), "w"):
                pass

        if dna_sequence:
            shutil.move(os.path.join(tmp, "amr_genes.fasta"), str(genes))
        else:
            with open(os.path.join(str(genes), "amr_genes.fasta"), "w"):
                pass

        if protein_sequence:
            shutil.move(os.path.join(tmp, "amr_proteins.fasta"), str(proteins))
        else:
            with open(os.path.join(str(proteins), "amr_proteins.fasta"), "w"):
                pass

    return annotations, mutations, genes, proteins
