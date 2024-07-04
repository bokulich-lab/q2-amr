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


def _path_maker(name):
    return str(name)


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

    def __init__(self, path, mode):
        super().__init__(path, mode)

        # Overwrite path maker methods for all file collections
        for var_name, var_value in vars(self.__class__).items():
            if isinstance(var_value, model.FileCollection):
                var_value.set_path_maker(_path_maker)
