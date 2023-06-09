# ----------------------------------------------------------------------------
# Copyright (c) 2022, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

from q2_types.feature_data import FeatureData, ProteinSequence, Sequence
from qiime2.core.type import Bool, Choices, Int, List, Range, Str
from qiime2.plugin import Citations, Plugin

from q2_amr import __version__
from q2_amr.card import annotate_card, fetch_card_db
from q2_amr.heatmap import heatmap
from q2_amr.types import (
    CARDAnnotation,
    CARDAnnotationDirectoryFormat,
    CARDAnnotationJSONFormat,
    CARDAnnotationTXTFormat,
    CARDDatabase,
    CARDDatabaseDirectoryFormat,
    CARDDatabaseFormat,
)

citations = Citations.load("citations.bib", package="q2_amr")

plugin = Plugin(
    name="amr",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-amr",
    package="q2_amr",
    description="This is a QIIME 2 plugin that annotates microbiome sequence data with "
    "antimicrobial resistance gene "
    "information from CARD.",
    short_description="This is a QIIME 2 plugin that annotates microbiome sequence data"
    " with antimicrobial resistance "
    "gene information from CARD.",
)
plugin.methods.register_function(
    function=fetch_card_db,
    inputs={},
    parameters={
        "version": Str
        % Choices(
            [
                "3.2.6",
                "3.2.5",
                "3.2.4",
                "3.2.3",
                "3.2.2",
                "3.2.1",
                "3.2.0",
                "3.1.4",
                "3.1.3",
                "3.1.2",
                "3.1.1",
                "3.1.0",
                "3.0.9",
                "3.0.8",
            ]
        )
    },
    outputs=[("card_db", CARDDatabase)],
    input_descriptions={},
    parameter_descriptions={
        "version": "Version of the CARD database to be downloaded."
    },
    output_descriptions={
        "card_db": "CARD database of resistance genes, their products and associated "
        "phenotypes."
    },
    name="Download CARD data.",
    description=("Downloads the CARD database from the CARD website."),
    citations=[citations["alcock_card_2023"]],
)

plugin.methods.register_function(
    function=annotate_card,
    inputs={"input_sequence": FeatureData[Sequence]},
    parameters={
        "alignment_tool": Str % Choices(["BLAST", "DIAMOND"]),
        "input_type": Str % Choices(["contig", "protein"]),
        "split_prodigal_jobs": Bool,
        "include_loose": Bool,
        "exclude_nudge": Bool,
        "low_quality": Bool,
        "num_threads": Int % Range(1, None),
    },
    outputs=[
        ("amr_annotations", CARDAnnotation),
        ("protein_annotation", FeatureData[ProteinSequence]),
        ("dna_annotation", FeatureData[Sequence]),
    ],
    input_descriptions={"input_sequence": "Sequences to be annotated with rgi."},
    parameter_descriptions={
        "alignment_tool": "Specify alignment tool BLAST or DIAMOND.",
        "input_type": "Specify data input type contig or protein.",
        "split_prodigal_jobs": "Run multiple prodigal jobs simultaneously for contigs "
        "in a fasta file.",
        "include_loose": "Include loose hits in addition to strict and perfect hits .",
        "exclude_nudge": "Include hits nudged from loose to strict hits.",
        "low_quality": "Use for short contigs to predict partial genes.",
        "num_threads": "Number of threads (CPUs) to use in the BLAST search.",
    },
    output_descriptions={
        "amr_annotations": "AMR Annotation as .txt and .json file.",
        "protein_annotation": "FASTA file with predicted protein sequences, ORF_ID and "
        "ARO accession in the Header.",
        "dna_annotation": "FASTA file with predicted dna sequences, ORF_ID and ARO "
        "accession in the Header.",
    },
    name="Annotation of sequence data with antimicrobial resistance gene information "
    "from CARD.",
    description=(
        "Annotation of sequence data with antimicrobial resistance gene information "
        "from CARD."
    ),
    citations=[citations["alcock_card_2023"]],
)


plugin.visualizers.register_function(
    function=heatmap,
    inputs={"amr_annotation": List[CARDAnnotation]},
    parameters={
        "cat": Str % Choices(["drug_class", "resistance_mechanism", "gene_family"]),
        "clus": Str % Choices(["samples", "genes", "both"]),
        "display": Str % Choices(["plain", "fill", "text"]),
        "frequency": Bool,
    },
    input_descriptions={"amr_annotation": "Sequences to be annotated with rgi."},
    parameter_descriptions={
        "cat": "The option to organize resistance genes based on a category.",
        "clus": "Specify data input type contig or protein.",
        "display": "Specify display options for categories",
        "frequency": "Represent samples based on resistance profile.",
    },
    name="Download CARD data.",
    description=("Downloads the CARD database from the CARD website."),
    citations=[citations["alcock_card_2023"]],
)

# Registrations
plugin.register_semantic_types(CARDDatabase, CARDAnnotation)

plugin.register_semantic_type_to_format(
    CARDDatabase, artifact_format=CARDDatabaseDirectoryFormat
)
plugin.register_semantic_type_to_format(
    CARDAnnotation, artifact_format=CARDAnnotationDirectoryFormat
)
plugin.register_formats(
    CARDAnnotationTXTFormat,
    CARDAnnotationDirectoryFormat,
    CARDAnnotationJSONFormat,
    CARDDatabaseFormat,
    CARDDatabaseDirectoryFormat,
)

importlib.import_module("q2_amr.types._transformer")
