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


class AMRFinderPlusDatabaseDirectoryFormat(model.DirectoryFormat):
    AMR_LIB = model.File("AMR.LIB", format=TextFormat)
    AMR_LIB_comp = model.FileCollection(r"AMR\.LIB\.h3.$", format=BinaryFormat)
    AMRProt = model.File("AMRProt", format=ProteinFASTAFormat)
    AMRProt_blast = model.FileCollection(r"AMRProt\.p..$", format=BinaryFormat)
    AMRProt_mutation = model.File("AMRProt-mutation.tab", format=TextFormat)
    AMRProt_suppress = model.File("AMRProt-suppress", format=TextFormat)
    AMRProt_susceptible = model.File("AMRProt-susceptible.tab", format=TextFormat)
    changes = model.File("changes.txt", format=TextFormat)
    db_version = model.File("database_format_version.txt", format=TextFormat)
    fam = model.File("fam.tab", format=TextFormat)
    taxgroup = model.File("taxgroup.tab", format=TextFormat)
    version = model.File("version.txt", format=TextFormat)
    AMR_DNA = model.FileCollection(
        r"^AMR_DNA-[a-zA-Z_]+$", format=MixedCaseDNAFASTAFormat
    )
    AMR_DNA_comp = model.FileCollection(
        r"^AMR_DNA-[a-zA-Z_]+\.n..$", format=BinaryFormat
    )
    AMR_CDS_comp = model.FileCollection(r"^AMR_CDS\.n..$", format=BinaryFormat)
    AMR_CDS = model.File("AMR_CDS", format=MixedCaseDNAFASTAFormat)
