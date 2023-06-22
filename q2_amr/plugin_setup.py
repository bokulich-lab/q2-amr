# ----------------------------------------------------------------------------
# Copyright (c) 2022, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

from q2_types.feature_table import FeatureTable, Frequency
from q2_types.per_sample_sequences import (
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData
from q2_types_genomics.per_sample_data import MAGs
from qiime2.core.type import Bool, Choices, Float, Int, Range, Str
from qiime2.plugin import Citations, Plugin

from q2_amr import __version__
from q2_amr.card.database import fetch_card_db
from q2_amr.card.heatmap import heatmap
from q2_amr.card.mags import annotate_mags_card
from q2_amr.card.reads import annotate_reads_card, visualize_annotation_stats
from q2_amr.types import (
    CARDAnnotationJSONFormat,
    CARDAnnotationTXTFormat,
    CARDDatabase,
    CARDDatabaseDirectoryFormat,
    CARDDatabaseFormat,
)
from q2_amr.types._format import (
    CARDAlleleAnnotationDirectoryFormat,
    CARDAlleleAnnotationFormat,
    CARDAnnotationDirectoryFormat,
    CARDAnnotationStatsFormat,
    CARDGeneAnnotationDirectoryFormat,
    CARDGeneAnnotationFormat,
)
from q2_amr.types._type import CARDAlleleAnnotation, CARDAnnotation, CARDGeneAnnotation

citations = Citations.load("citations.bib", package="q2_amr")

plugin = Plugin(
    name="amr",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-amr",
    package="q2_amr",
    description="This is a QIIME 2 plugin that annotates microbiome sequence data with "
    "antimicrobial resistance gene information from CARD.",
    short_description="This is a QIIME 2 plugin that annotates microbiome sequence "
    "data with antimicrobial resistance gene information from CARD.",
)
plugin.methods.register_function(
    function=fetch_card_db,
    inputs={},
    parameters={},
    outputs=[("card_db", CARDDatabase)],
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={
        "card_db": "CARD database of resistance genes, their products and associated "
        "phenotypes."
    },
    name="Download CARD data.",
    description=("Downloads the CARD database from the CARD website."),
    citations=[citations["alcock_card_2023"]],
)

plugin.methods.register_function(
    function=annotate_mags_card,
    inputs={"mag": SampleData[MAGs], "card_db": CARDDatabase},
    parameters={
        "alignment_tool": Str % Choices(["BLAST", "DIAMOND"]),
        "input_type": Str % Choices(["contig", "protein"]),
        "split_prodigal_jobs": Bool,
        "include_loose": Bool,
        "include_nudge": Bool,
        "low_quality": Bool,
        "num_threads": Int % Range(1, 9),
    },
    outputs=[("amr_annotations", SampleData[CARDAnnotation])],
    input_descriptions={
        "mag": "MAG to be annotated with CARD.",
        "card_db": "CARD Database.",
    },
    parameter_descriptions={
        "alignment_tool": "Specify alignment tool BLAST or DIAMOND.",
        "input_type": "Specify data input type contig or protein.",
        "split_prodigal_jobs": "Run multiple prodigal jobs simultaneously for contigs "
        "in a fasta file.",
        "include_loose": "Include loose hits in addition to strict and perfect hits .",
        "include_nudge": "Include hits nudged from loose to strict hits.",
        "low_quality": "Use for short contigs to predict partial genes.",
        "num_threads": "Number of threads (CPUs) to use in the BLAST search.",
    },
    output_descriptions={"amr_annotations": "AMR Annotation as .txt and .json file."},
    name="Annotate MAGs with antimicrobial resistance gene information from CARD.",
    description="Annotate MAGs with antimicrobial resistance gene information from "
    "CARD.",
    citations=[citations["alcock_card_2023"]],
)


plugin.methods.register_function(
    function=annotate_reads_card,
    inputs={
        "reads": SampleData[PairedEndSequencesWithQuality | SequencesWithQuality],
        "card_db": CARDDatabase,
    },
    parameters={
        "aligner": Str % Choices(["kma", "bowtie2", "bwa"]),
        "include_baits": Bool,
        "mapq": Float % Range(0, None, inclusive_start=True),
        "mapped": Float % Range(0, None, inclusive_start=True),
        "coverage": Float % Range(0, None, inclusive_start=True),
        "threads": Int % Range(0, None, inclusive_start=False),
    },
    outputs=[
        ("amr_allele_annotation", CARDAlleleAnnotation),
        ("amr_gene_annotation", CARDGeneAnnotation),
        ("allele_feature_table", FeatureTable[Frequency]),
        ("gene_feature_table", FeatureTable[Frequency]),
    ],
    input_descriptions={
        "reads": "Paired or single end metagenomic reads.",
        "card_db": "CARD Database",
    },
    parameter_descriptions={
        "aligner": "Specify alignment tool.",
        "include_baits": "Include baits.",
        "mapq": "Filter reads based on MAPQ score.",
        "mapped": "Filter reads based on mapped reads.",
        "coverage": "Filter reads based on coverage of reference sequence.",
        "threads": "Number of threads (CPUs) to use.",
    },
    output_descriptions={
        "amr_allele_annotation": "AMR annotation mapped on alleles.",
        "amr_gene_annotation": "AMR annotation mapped on genes.",
        "allele_feature_table": "Samples combined into one frequency count table.",
        "gene_feature_table": "Samples combined into one frequency count table.",
    },
    name="Annotate metagenomic reads with antimicrobial resistance gene information "
    "from CARD.",
    description="Annotate metagenomic reads with antimicrobial resistance gene "
    "information from CARD.",
    citations=[citations["alcock_card_2023"]],
)

plugin.visualizers.register_function(
    function=heatmap,
    inputs={"amr_annotation": CARDAnnotation},
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
    name="Create heatmap from annotate_mags_card output.",
    description=("Create heatmap from annotate_mags_card output."),
    citations=[citations["alcock_card_2023"]],
)

plugin.visualizers.register_function(
    function=visualize_annotation_stats,
    inputs={"amr_reads_annotation": CARDGeneAnnotation | CARDAlleleAnnotation},
    parameters={},
    input_descriptions={
        "amr_reads_annotation": "AMR annotation mapped on alleles or genes."
    },
    parameter_descriptions={},
    name="Visualize mapping statistics.",
    description="Visualize mapping statistics.",
    citations=[citations["alcock_card_2023"]],
)

# Registrations
plugin.register_semantic_types(
    CARDDatabase, CARDAnnotation, CARDAlleleAnnotation, CARDGeneAnnotation
)

plugin.register_semantic_type_to_format(
    CARDDatabase, artifact_format=CARDDatabaseDirectoryFormat
)
plugin.register_semantic_type_to_format(
    CARDAnnotation, artifact_format=CARDAnnotationDirectoryFormat
)
plugin.register_semantic_type_to_format(
    CARDAlleleAnnotation, artifact_format=CARDAlleleAnnotationDirectoryFormat
)
plugin.register_semantic_type_to_format(
    CARDGeneAnnotation, artifact_format=CARDGeneAnnotationDirectoryFormat
)

plugin.register_formats(
    CARDAnnotationTXTFormat,
    CARDAnnotationJSONFormat,
    CARDAnnotationDirectoryFormat,
    CARDDatabaseFormat,
    CARDDatabaseDirectoryFormat,
    CARDAlleleAnnotationFormat,
    CARDGeneAnnotationFormat,
    CARDAnnotationStatsFormat,
    CARDAlleleAnnotationDirectoryFormat,
    CARDGeneAnnotationDirectoryFormat,
)

importlib.import_module("q2_amr.types._transformer")
