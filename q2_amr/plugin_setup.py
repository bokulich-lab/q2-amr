# ----------------------------------------------------------------------------
# Copyright (c) 2022, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

from q2_types.per_sample_sequences import PairedEndSequencesWithQuality
from q2_types.sample_data import SampleData

from q2_amr.types import CARDDatabase, CARDDatabaseDirectoryFormat, CARDAnnotationTXTFormat, \
    CARDDatabaseFormat, CARDAnnotationJSONFormat
from q2_types.feature_data import Sequence, FeatureData, ProteinSequence
from qiime2.core.type import Str, Choices, Bool, Int, Range, Float

from q2_amr.card import fetch_card_db, annotate_card, heatmap, bwt  # heatmap
from qiime2.plugin import Citations, Plugin

from q2_amr import __version__
from q2_amr.types._format import CARDAnnotationDirectoryFormat
from q2_amr.types._type import CARDAnnotation

citations = Citations.load("citations.bib", package="q2_amr")

plugin = Plugin(
    name="amr",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-amr",
    package="q2_amr",
    description="This is a QIIME 2 plugin that annotates microbiome sequence data with antimicrobial resistance gene "
                "information from CARD.",
    short_description="This is a QIIME 2 plugin that annotates microbiome sequence data with antimicrobial resistance "
                      "gene information from CARD.",
)
plugin.methods.register_function(
    function=fetch_card_db,
    inputs={},
    parameters={},
    outputs=[('card_db', CARDDatabase)],
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={
        'card_db': 'CARD database of resistance genes, their products and associated '
                   'phenotypes.'},
    name='Download CARD data.',
    description=('Downloads the CARD database from the CARD website.'),
    citations=[citations['alcock_card_2023']]
)

plugin.methods.register_function(
    function=annotate_card,
    inputs={'input_sequence': FeatureData[Sequence]},
    parameters={'alignment_tool': Str % Choices(['BLAST', 'DIAMOND']),
                'input_type': Str % Choices(['contig', 'protein']),
                'split_prodigal_jobs': Bool,
                'include_loose': Bool,
                'exclude_nudge': Bool,
                'low_quality': Bool,
                'num_threads': Int % Range(1, None)},
    outputs=[('amr_annotations', CARDAnnotation),
             ('protein_annotation', FeatureData[ProteinSequence]),
             ('dna_annotation', FeatureData[Sequence])],
    input_descriptions={'input_sequence': 'Sequences to be annotated with rgi.'},
    parameter_descriptions={
        'alignment_tool': 'Specify alignment tool BLAST or DIAMOND.',
        'input_type': 'Specify data input type contig or protein.',
        'split_prodigal_jobs': 'Run multiple prodigal jobs simultaneously for contigs in a fasta file.',
        'include_loose': 'Include loose hits in addition to strict and perfect hits .',
        'exclude_nudge': 'Include hits nudged from loose to strict hits.',
        'low_quality': 'Use for short contigs to predict partial genes.',
        'num_threads': 'Number of threads (CPUs) to use in the BLAST search.'},
    output_descriptions={
        'amr_annotations': 'AMR Annotation as .txt and .json file.',
        'protein_annotation': 'FASTA file with predicted protein sequences, ORF_ID and ARO accession in the Header.',
        'dna_annotation': 'FASTA file with predicted dna sequences, ORF_ID and ARO accession in the Header.'},
    name='Annotation of sequence data with antimicrobial resistance gene information from CARD.',
    description=('Annotation of sequence data with antimicrobial resistance gene information from CARD.'),
    citations=[citations['alcock_card_2023']]
)

plugin.visualizers.register_function(
    function=heatmap,
    inputs={'amr_annotation_json': CARDAnnotation},
    parameters={},
    input_descriptions={'amr_annotation_json': 'Sequences to be annotated with rgi.'},
    parameter_descriptions={},
    name='Download CARD data.',
    description=('Downloads the CARD database from the CARD website.'),
    citations=[citations['alcock_card_2023']]
)

plugin.visualizers.register_function(
    function=bwt,
    inputs={'reads': SampleData[PairedEndSequencesWithQuality],
            'card_db': CARDDatabase},
    parameters={'aligner': Str % Choices(['kma', 'bowtie2', 'bwa']),
                'local': Bool,
                'include_baits': Bool,
                'mapq': Int % Range(1, None),
                'mapped': Int % Range(1, None),
                'coverage': Float % Range(0, None, inclusive_start=False),
                'threads': Int},
    input_descriptions={'reads': 'Paired end metagenomic reads.',
                        'card_db': 'CARD Database'},
    parameter_descriptions={
        'aligner': 'Specify alignment tool kma, bowtie2 or bwa.',
        'local': 'Use local database (default: uses database in executable directory).',
        'include_baits': 'Include baits.',
        'mapq': 'Filter reads based on MAPQ score.',
        'mapped': 'Filter reads based on mapped reads.',
        'coverage': 'Filter reads based on coverage of reference sequence.',
        'threads': 'Number of threads (CPUs) to use.'},
    name='Annotation of metagenomic reads with antimicrobial resistance gene information from CARD.',
    description=('Annotation of metagenomic reads with antimicrobial resistance gene information from CARD.'),
    citations=[citations['alcock_card_2023']]
)

#Registrations
plugin.register_semantic_types(CARDDatabase, CARDAnnotation)

plugin.register_semantic_type_to_format(
    CARDDatabase,
    artifact_format=CARDDatabaseDirectoryFormat)
plugin.register_semantic_type_to_format(
    CARDAnnotation,
    artifact_format=CARDAnnotationDirectoryFormat)
plugin.register_formats(CARDAnnotationTXTFormat, CARDAnnotationDirectoryFormat,
                        CARDAnnotationJSONFormat,
                        CARDDatabaseFormat, CARDDatabaseDirectoryFormat)

importlib.import_module('q2_amr.types._transformer')
