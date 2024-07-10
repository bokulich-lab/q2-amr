# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from q2_types.feature_data import MixedCaseDNAFASTAFormat, ProteinFASTAFormat
from q2_types.per_sample_sequences._format import MultiDirValidationMixin
from qiime2.core.exceptions import ValidationError
from qiime2.plugin import model


class TextFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


class BinaryFormat(model.BinaryFileFormat):
    def _validate_(self, level):
        pass


class AMRFinderPlusDatabaseDirFmt(model.DirectoryFormat):
    amr_lib = model.File("AMR.LIB", format=TextFormat)
    amr_lib_comp = model.FileCollection(r"^AMR\.LIB\.h3.$", format=BinaryFormat)
    amrprot = model.File("AMRProt", format=ProteinFASTAFormat)
    amrprot_blast = model.FileCollection(r"^AMRProt\.p..$", format=BinaryFormat)
    amrprot_mutation = model.File("AMRProt-mutation.tab", format=TextFormat)
    amrprot_suppress = model.File("AMRProt-suppress", format=TextFormat)
    amrprot_susceptible = model.File("AMRProt-susceptible.tab", format=TextFormat)
    changes = model.File("changes.txt", format=TextFormat)
    db_version = model.File("database_format_version.txt", format=TextFormat)
    fam = model.File("fam.tab", format=TextFormat)
    taxgroup = model.File("taxgroup.tab", format=TextFormat)
    version = model.File("version.txt", format=TextFormat)
    amr_dna = model.FileCollection(
        r"^AMR_DNA-[a-zA-Z_]+$", format=MixedCaseDNAFASTAFormat
    )
    amr_dna_comp = model.FileCollection(
        r"^AMR_DNA-[a-zA-Z_]+\.n..$", format=BinaryFormat
    )
    amr_dna_tab = model.FileCollection(r"^AMR_DNA-[a-zA-Z_]+\.tab$", format=TextFormat)
    amr_cds_comp = model.FileCollection(r"^AMR_CDS\.n..$", format=BinaryFormat)
    amr_cds = model.File("AMR_CDS", format=MixedCaseDNAFASTAFormat)

    @amr_lib_comp.set_path_maker
    def amr_lib_comp_path_maker(self):
        return r"^AMR\.LIB\.h3.$"

    @amrprot_blast.set_path_maker
    def amrprot_blast_path_maker(self):
        return r"^AMRProt\.p..$"

    @amr_dna.set_path_maker
    def amr_dna_path_maker(self):
        return r"^AMR_DNA-[a-zA-Z_]+$"

    @amr_dna_comp.set_path_maker
    def amr_dna_comp_path_maker(self):
        return r"^AMR_DNA-[a-zA-Z_]+\.n..$"

    @amr_cds_comp.set_path_maker
    def amr_cds_comp_path_maker(self):
        return r"^AMR_CDS\.n..$"

    @amr_dna_tab.set_path_maker
    def amr_dna_tab_path_maker(self):
        return r"^AMR_DNA-[a-zA-Z_]+\.tab$"


class ARMFinderPlusAnnotationFormat(model.TextFileFormat):
    def _validate(self, n_records=None):
        header_coordinates = [
            "Protein identifier",
            "Contig id",
            "Start",
            "Stop",
            "Strand",
            "Gene symbol",
            "Sequence name",
            "Scope",
            "Element type",
            "Element subtype",
            "Class",
            "Subclass",
            "Method",
            "Target length",
            "Reference sequence length",
            "% Coverage of reference sequence",
            "% Identity to reference sequence",
            "Alignment length",
            "Accession of closest sequence",
            "Name of closest sequence",
            "HMM id",
            "HMM description",
            "Hierarchy node",
        ]
        header = header_coordinates[:1] + header_coordinates[5:]
        header_obs = pd.read_csv(str(self), sep="\t", nrows=0).columns.tolist()
        if header != header_obs and header_coordinates != header_obs:
            raise ValidationError(
                "Header line does not match ARMFinderPlusAnnotation format. Must "
                "consist of the following values: "
                + ", ".join(header_coordinates)
                + ".\nWhile Contig id, Start, Stop and Strand are optional."
                + ".\n\nFound instead: "
                + ", ".join(header_obs)
            )

    def _validate_(self, level):
        self._validate()


class ARMFinderPlusAnnotationsDirFmt(MultiDirValidationMixin, model.DirectoryFormat):
    annotation = model.FileCollection(
        r".+amr_(annotations|mutations)\.tsv$", format=ARMFinderPlusAnnotationFormat
    )

    @annotation.set_path_maker
    def annotation_path_maker(self, sample_id, mag_id):
        prefix = f"{sample_id}/{mag_id}_" if mag_id else f"{sample_id}/"
        return f"{prefix}amr_annotations.tsv"


ARMFinderPlusAnnotationDirFmt = model.SingleFileDirectoryFormat(
    "ARMFinderPlusAnnotationDirFmt",
    r"amr_(annotations|mutations)\.tsv$",
    ARMFinderPlusAnnotationFormat,
)
