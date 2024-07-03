# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from q2_types.feature_data import MixedCaseDNAFASTAFormat, ProteinFASTAFormat
from qiime2.core.exceptions import ValidationError
from qiime2.plugin import model


class TextFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


class BinaryFormat(model.BinaryFileFormat):
    def _validate_(self, level):
        pass


class AMRFinderPlusDatabaseDirectoryFormat(model.DirectoryFormat):
    AMR_LIB = model.File("AMR.LIB", format=TextFormat)
    AMRProt = model.File("AMRProt", format=ProteinFASTAFormat)
    AMRProt_mutation = model.File("AMRProt-mutation.tab", format=TextFormat)
    AMRProt_suppress = model.File("AMRProt-suppress", format=TextFormat)
    AMRProt_susceptible = model.File("AMRProt-susceptible.tab", format=TextFormat)
    changes = model.File("changes.txt", format=TextFormat)
    db_version = model.File("database_format_version.txt", format=TextFormat)
    fam = model.File("fam.tab", format=TextFormat)
    taxgroup = model.File("taxgroup.tab", format=TextFormat)
    version = model.File("version.txt", format=TextFormat)
    AMR_CDS = model.File("AMR_CDS", format=MixedCaseDNAFASTAFormat)
    AMR_LIB_h3f = model.File("AMR.LIB.h3f", format=BinaryFormat)
    AMR_LIB_h3i = model.File("AMR.LIB.h3i", format=BinaryFormat)
    AMR_LIB_h3m = model.File("AMR.LIB.h3m", format=BinaryFormat)
    AMR_LIB_h3p = model.File("AMR.LIB.h3p", format=BinaryFormat)
    AMRProt_pdb = model.File("AMRProt.pdb", format=BinaryFormat)
    AMRProt_phr = model.File("AMRProt.phr", format=BinaryFormat)
    AMRProt_pin = model.File("AMRProt.pin", format=BinaryFormat)
    AMRProt_pjs = model.File("AMRProt.pjs", format=BinaryFormat)
    AMRProt_pot = model.File("AMRProt.pot", format=BinaryFormat)
    AMRProt_psq = model.File("AMRProt.psq", format=BinaryFormat)
    AMRProt_ptf = model.File("AMRProt.ptf", format=BinaryFormat)
    AMRProt_pto = model.File("AMRProt.pto", format=BinaryFormat)
    AMR_CDS_ndb = model.File("AMR_CDS.ndb", format=BinaryFormat)
    AMR_CDS_nhr = model.File("AMR_CDS.nhr", format=BinaryFormat)
    AMR_CDS_nin = model.File("AMR_CDS.nin", format=BinaryFormat)
    AMR_CDS_njs = model.File("AMR_CDS.njs", format=BinaryFormat)
    AMR_CDS_not = model.File("AMR_CDS.not", format=BinaryFormat)
    AMR_CDS_nsq = model.File("AMR_CDS.nsq", format=BinaryFormat)
    AMR_CDS_ntf = model.File("AMR_CDS.ntf", format=BinaryFormat)
    AMR_CDS_nto = model.File("AMR_CDS.nto", format=BinaryFormat)
    AMR_DNA_Acinetobacter_baumannii_ndb = model.File(
        "AMR_DNA-Acinetobacter_baumannii.ndb", format=BinaryFormat
    )
    AMR_DNA_Acinetobacter_baumannii_nhr = model.File(
        "AMR_DNA-Acinetobacter_baumannii.nhr", format=BinaryFormat
    )
    AMR_DNA_Acinetobacter_baumannii_nin = model.File(
        "AMR_DNA-Acinetobacter_baumannii.nin", format=BinaryFormat
    )
    AMR_DNA_Acinetobacter_baumannii_njs = model.File(
        "AMR_DNA-Acinetobacter_baumannii.njs", format=BinaryFormat
    )
    AMR_DNA_Acinetobacter_baumannii_not = model.File(
        "AMR_DNA-Acinetobacter_baumannii.not", format=BinaryFormat
    )
    AMR_DNA_Acinetobacter_baumannii_nsq = model.File(
        "AMR_DNA-Acinetobacter_baumannii.nsq", format=BinaryFormat
    )
    AMR_DNA_Acinetobacter_baumannii_ntf = model.File(
        "AMR_DNA-Acinetobacter_baumannii.ntf", format=BinaryFormat
    )
    AMR_DNA_Acinetobacter_baumannii_nto = model.File(
        "AMR_DNA-Acinetobacter_baumannii.nto", format=BinaryFormat
    )
    AMR_DNA_Campylobacter_ndb = model.File(
        "AMR_DNA-Campylobacter.ndb", format=BinaryFormat
    )
    AMR_DNA_Campylobacter_nhr = model.File(
        "AMR_DNA-Campylobacter.nhr", format=BinaryFormat
    )
    AMR_DNA_Campylobacter_nin = model.File(
        "AMR_DNA-Campylobacter.nin", format=BinaryFormat
    )
    AMR_DNA_Campylobacter_njs = model.File(
        "AMR_DNA-Campylobacter.njs", format=BinaryFormat
    )
    AMR_DNA_Campylobacter_not = model.File(
        "AMR_DNA-Campylobacter.not", format=BinaryFormat
    )
    AMR_DNA_Campylobacter_nsq = model.File(
        "AMR_DNA-Campylobacter.nsq", format=BinaryFormat
    )
    AMR_DNA_Campylobacter_ntf = model.File(
        "AMR_DNA-Campylobacter.ntf", format=BinaryFormat
    )
    AMR_DNA_Campylobacter_nto = model.File(
        "AMR_DNA-Campylobacter.nto", format=BinaryFormat
    )
    AMR_DNA_Clostridioides_difficile_ndb = model.File(
        "AMR_DNA-Clostridioides_difficile.ndb", format=BinaryFormat
    )
    AMR_DNA_Clostridioides_difficile_nhr = model.File(
        "AMR_DNA-Clostridioides_difficile.nhr", format=BinaryFormat
    )
    AMR_DNA_Clostridioides_difficile_nin = model.File(
        "AMR_DNA-Clostridioides_difficile.nin", format=BinaryFormat
    )
    AMR_DNA_Clostridioides_difficile_njs = model.File(
        "AMR_DNA-Clostridioides_difficile.njs", format=BinaryFormat
    )
    AMR_DNA_Clostridioides_difficile_not = model.File(
        "AMR_DNA-Clostridioides_difficile.not", format=BinaryFormat
    )
    AMR_DNA_Clostridioides_difficile_nsq = model.File(
        "AMR_DNA-Clostridioides_difficile.nsq", format=BinaryFormat
    )
    AMR_DNA_Clostridioides_difficile_ntf = model.File(
        "AMR_DNA-Clostridioides_difficile.ntf", format=BinaryFormat
    )
    AMR_DNA_Clostridioides_difficile_nto = model.File(
        "AMR_DNA-Clostridioides_difficile.nto", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecalis_ndb = model.File(
        "AMR_DNA-Enterococcus_faecalis.ndb", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecalis_nhr = model.File(
        "AMR_DNA-Enterococcus_faecalis.nhr", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecalis_nin = model.File(
        "AMR_DNA-Enterococcus_faecalis.nin", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecalis_njs = model.File(
        "AMR_DNA-Enterococcus_faecalis.njs", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecalis_not = model.File(
        "AMR_DNA-Enterococcus_faecalis.not", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecalis_nsq = model.File(
        "AMR_DNA-Enterococcus_faecalis.nsq", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecalis_ntf = model.File(
        "AMR_DNA-Enterococcus_faecalis.ntf", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecalis_nto = model.File(
        "AMR_DNA-Enterococcus_faecalis.nto", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecium_ndb = model.File(
        "AMR_DNA-Enterococcus_faecium.ndb", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecium_nhr = model.File(
        "AMR_DNA-Enterococcus_faecium.nhr", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecium_nin = model.File(
        "AMR_DNA-Enterococcus_faecium.nin", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecium_njs = model.File(
        "AMR_DNA-Enterococcus_faecium.njs", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecium_not = model.File(
        "AMR_DNA-Enterococcus_faecium.not", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecium_nsq = model.File(
        "AMR_DNA-Enterococcus_faecium.nsq", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecium_ntf = model.File(
        "AMR_DNA-Enterococcus_faecium.ntf", format=BinaryFormat
    )
    AMR_DNA_Enterococcus_faecium_nto = model.File(
        "AMR_DNA-Enterococcus_faecium.nto", format=BinaryFormat
    )
    AMR_DNA_Escherichia_ndb = model.File("AMR_DNA-Escherichia.ndb", format=BinaryFormat)
    AMR_DNA_Escherichia_nhr = model.File("AMR_DNA-Escherichia.nhr", format=BinaryFormat)
    AMR_DNA_Escherichia_nin = model.File("AMR_DNA-Escherichia.nin", format=BinaryFormat)
    AMR_DNA_Escherichia_njs = model.File("AMR_DNA-Escherichia.njs", format=BinaryFormat)
    AMR_DNA_Escherichia_not = model.File("AMR_DNA-Escherichia.not", format=BinaryFormat)
    AMR_DNA_Escherichia_nsq = model.File("AMR_DNA-Escherichia.nsq", format=BinaryFormat)
    AMR_DNA_Escherichia_ntf = model.File("AMR_DNA-Escherichia.ntf", format=BinaryFormat)
    AMR_DNA_Escherichia_nto = model.File("AMR_DNA-Escherichia.nto", format=BinaryFormat)
    AMR_DNA_Klebsiella_oxytoca_ndb = model.File(
        "AMR_DNA-Klebsiella_oxytoca.ndb", format=BinaryFormat
    )
    AMR_DNA_Klebsiella_oxytoca_nhr = model.File(
        "AMR_DNA-Klebsiella_oxytoca.nhr", format=BinaryFormat
    )
    AMR_DNA_Klebsiella_oxytoca_nin = model.File(
        "AMR_DNA-Klebsiella_oxytoca.nin", format=BinaryFormat
    )
    AMR_DNA_Klebsiella_oxytoca_njs = model.File(
        "AMR_DNA-Klebsiella_oxytoca.njs", format=BinaryFormat
    )
    AMR_DNA_Klebsiella_oxytoca_not = model.File(
        "AMR_DNA-Klebsiella_oxytoca.not", format=BinaryFormat
    )
    AMR_DNA_Klebsiella_oxytoca_nsq = model.File(
        "AMR_DNA-Klebsiella_oxytoca.nsq", format=BinaryFormat
    )
    AMR_DNA_Klebsiella_oxytoca_ntf = model.File(
        "AMR_DNA-Klebsiella_oxytoca.ntf", format=BinaryFormat
    )
    AMR_DNA_Klebsiella_oxytoca_nto = model.File(
        "AMR_DNA-Klebsiella_oxytoca.nto", format=BinaryFormat
    )
    AMR_DNA_Neisseria_gonorrhoeae_ndb = model.File(
        "AMR_DNA-Neisseria_gonorrhoeae.ndb", format=BinaryFormat
    )
    AMR_DNA_Neisseria_gonorrhoeae_nhr = model.File(
        "AMR_DNA-Neisseria_gonorrhoeae.nhr", format=BinaryFormat
    )
    AMR_DNA_Neisseria_gonorrhoeae_nin = model.File(
        "AMR_DNA-Neisseria_gonorrhoeae.nin", format=BinaryFormat
    )
    AMR_DNA_Neisseria_gonorrhoeae_njs = model.File(
        "AMR_DNA-Neisseria_gonorrhoeae.njs", format=BinaryFormat
    )
    AMR_DNA_Neisseria_gonorrhoeae_not = model.File(
        "AMR_DNA-Neisseria_gonorrhoeae.not", format=BinaryFormat
    )
    AMR_DNA_Neisseria_gonorrhoeae_nsq = model.File(
        "AMR_DNA-Neisseria_gonorrhoeae.nsq", format=BinaryFormat
    )
    AMR_DNA_Neisseria_gonorrhoeae_ntf = model.File(
        "AMR_DNA-Neisseria_gonorrhoeae.ntf", format=BinaryFormat
    )
    AMR_DNA_Neisseria_gonorrhoeae_nto = model.File(
        "AMR_DNA-Neisseria_gonorrhoeae.nto", format=BinaryFormat
    )
    AMR_DNA_Salmonella_nhr = model.File("AMR_DNA-Salmonella.nhr", format=BinaryFormat)
    AMR_DNA_Salmonella_ndb = model.File("AMR_DNA-Salmonella.ndb", format=BinaryFormat)
    AMR_DNA_Salmonella_nin = model.File("AMR_DNA-Salmonella.nin", format=BinaryFormat)
    AMR_DNA_Salmonella_njs = model.File("AMR_DNA-Salmonella.njs", format=BinaryFormat)
    AMR_DNA_Salmonella_not = model.File("AMR_DNA-Salmonella.not", format=BinaryFormat)
    AMR_DNA_Salmonella_nsq = model.File("AMR_DNA-Salmonella.nsq", format=BinaryFormat)
    AMR_DNA_Salmonella_ntf = model.File("AMR_DNA-Salmonella.ntf", format=BinaryFormat)
    AMR_DNA_Salmonella_nto = model.File("AMR_DNA-Salmonella.nto", format=BinaryFormat)
    AMR_DNA_Staphylococcus_aureus_ndb = model.File(
        "AMR_DNA-Staphylococcus_aureus.ndb", format=BinaryFormat
    )
    AMR_DNA_Staphylococcus_aureus_nhr = model.File(
        "AMR_DNA-Staphylococcus_aureus.nhr", format=BinaryFormat
    )
    AMR_DNA_Staphylococcus_aureus_nin = model.File(
        "AMR_DNA-Staphylococcus_aureus.nin", format=BinaryFormat
    )
    AMR_DNA_Staphylococcus_aureus_njs = model.File(
        "AMR_DNA-Staphylococcus_aureus.njs", format=BinaryFormat
    )
    AMR_DNA_Staphylococcus_aureus_not = model.File(
        "AMR_DNA-Staphylococcus_aureus.not", format=BinaryFormat
    )
    AMR_DNA_Staphylococcus_aureus_nsq = model.File(
        "AMR_DNA-Staphylococcus_aureus.nsq", format=BinaryFormat
    )
    AMR_DNA_Staphylococcus_aureus_ntf = model.File(
        "AMR_DNA-Staphylococcus_aureus.ntf", format=BinaryFormat
    )
    AMR_DNA_Staphylococcus_aureus_nto = model.File(
        "AMR_DNA-Staphylococcus_aureus.nto", format=BinaryFormat
    )
    AMR_DNA_Streptococcus_pneumoniae_ndb = model.File(
        "AMR_DNA-Streptococcus_pneumoniae.ndb", format=BinaryFormat
    )
    AMR_DNA_Streptococcus_pneumoniae_nhr = model.File(
        "AMR_DNA-Streptococcus_pneumoniae.nhr", format=BinaryFormat
    )
    AMR_DNA_Streptococcus_pneumoniae_nin = model.File(
        "AMR_DNA-Streptococcus_pneumoniae.nin", format=BinaryFormat
    )
    AMR_DNA_Streptococcus_pneumoniae_njs = model.File(
        "AMR_DNA-Streptococcus_pneumoniae.njs", format=BinaryFormat
    )
    AMR_DNA_Streptococcus_pneumoniae_not = model.File(
        "AMR_DNA-Streptococcus_pneumoniae.not", format=BinaryFormat
    )
    AMR_DNA_Streptococcus_pneumoniae_nsq = model.File(
        "AMR_DNA-Streptococcus_pneumoniae.nsq", format=BinaryFormat
    )
    AMR_DNA_Streptococcus_pneumoniae_ntf = model.File(
        "AMR_DNA-Streptococcus_pneumoniae.ntf", format=BinaryFormat
    )
    AMR_DNA_Streptococcus_pneumoniae_nto = model.File(
        "AMR_DNA-Streptococcus_pneumoniae.nto", format=BinaryFormat
    )
    AMR_DNA_Acinetobacter_baumannii = model.File(
        "AMR_DNA-Acinetobacter_baumannii", format=MixedCaseDNAFASTAFormat
    )
    AMR_DNA_Campylobacter = model.File(
        "AMR_DNA-Campylobacter", format=MixedCaseDNAFASTAFormat
    )
    AMR_DNA_Clostridioides_difficile = model.File(
        "AMR_DNA-Clostridioides_difficile", format=MixedCaseDNAFASTAFormat
    )
    AMR_DNA_Enterococcus_faecalis = model.File(
        "AMR_DNA-Enterococcus_faecalis", format=MixedCaseDNAFASTAFormat
    )
    AMR_DNA_Enterococcus_faecium = model.File(
        "AMR_DNA-Enterococcus_faecium", format=MixedCaseDNAFASTAFormat
    )
    AMR_DNA_Escherichia = model.File(
        "AMR_DNA-Escherichia", format=MixedCaseDNAFASTAFormat
    )
    AMR_DNA_Klebsiella_oxytoca = model.File(
        "AMR_DNA-Klebsiella_oxytoca", format=MixedCaseDNAFASTAFormat
    )
    AMR_DNA_Neisseria_gonorrhoeae = model.File(
        "AMR_DNA-Neisseria_gonorrhoeae", format=MixedCaseDNAFASTAFormat
    )
    AMR_DNA_Salmonella = model.File(
        "AMR_DNA-Salmonella", format=MixedCaseDNAFASTAFormat
    )
    AMR_DNA_Staphylococcus_aureus = model.File(
        "AMR_DNA-Staphylococcus_aureus", format=MixedCaseDNAFASTAFormat
    )
    AMR_DNA_Streptococcus_pneumoniae = model.File(
        "AMR_DNA-Streptococcus_pneumoniae", format=MixedCaseDNAFASTAFormat
    )
    AMR_DNA_Acinetobacter_baumannii_tab = model.File(
        "AMR_DNA-Acinetobacter_baumannii.tab", format=TextFormat
    )
    AMR_DNA_Campylobacter_tab = model.File(
        "AMR_DNA-Campylobacter.tab", format=TextFormat
    )
    AMR_DNA_Clostridioides_difficile_tab = model.File(
        "AMR_DNA-Clostridioides_difficile.tab", format=TextFormat
    )
    AMR_DNA_Enterococcus_faecalis_tab = model.File(
        "AMR_DNA-Enterococcus_faecalis.tab", format=TextFormat
    )
    AMR_DNA_Enterococcus_faecium_tab = model.File(
        "AMR_DNA-Enterococcus_faecium.tab", format=TextFormat
    )
    AMR_DNA_Klebsiella_oxytoca_tab = model.File(
        "AMR_DNA-Klebsiella_oxytoca.tab", format=TextFormat
    )
    AMR_DNA_Neisseria_gonorrhoeae_tab = model.File(
        "AMR_DNA-Neisseria_gonorrhoeae.tab", format=TextFormat
    )
    AMR_DNA_Salmonella_tab = model.File("AMR_DNA-Salmonella.tab", format=TextFormat)
    AMR_DNA_Staphylococcus_aureus_tab = model.File(
        "AMR_DNA-Staphylococcus_aureus.tab", format=TextFormat
    )
    AMR_DNA_Streptococcus_pneumoniae_tab = model.File(
        "AMR_DNA-Streptococcus_pneumoniae.tab", format=TextFormat
    )
    AMR_DNA_Escherichia_tab = model.File("AMR_DNA-Escherichia.tab", format=TextFormat)


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
