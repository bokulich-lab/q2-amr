# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_types.feature_data import MixedCaseDNAFASTAFormat, ProteinFASTAFormat
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
