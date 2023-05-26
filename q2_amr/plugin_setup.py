# ----------------------------------------------------------------------------
# Copyright (c) 2022, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

from q2_types.sample_data import SampleData
from q2_types_genomics.per_sample_data import MAGs

from q2_amr.types import CARDDatabase, CARDDatabaseDirectoryFormat, CARDAnnotationTXTFormat, \
    CARDDatabaseFormat, CARDAnnotationJSONFormat
from q2_types.feature_data import Sequence, FeatureData, ProteinSequence
from qiime2.core.type import Str, Choices, Bool, Int, Range

from q2_amr.card import fetch_card_db, annotate_mags_card, heatmap  # heatmap
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
    function=annotate_mags_card,
    inputs={'mag': SampleData[MAGs],
            'card_db': CARDDatabase},
    parameters={'alignment_tool': Str % Choices(['BLAST', 'DIAMOND']),
                'input_type': Str % Choices(['contig', 'protein']),
                'split_prodigal_jobs': Bool,
                'include_loose': Bool,
                'include_nudge': Bool,
                'low_quality': Bool,
                'num_threads': Int % Range(1, None)},
    outputs=[('amr_annotations', SampleData[CARDAnnotation])],
    input_descriptions={'mag': 'MAG to be annotated with CARD.',
                        'card_db': 'CARD Database.'},
    parameter_descriptions={
        'alignment_tool': 'Specify alignment tool BLAST or DIAMOND.',
        'input_type': 'Specify data input type contig or protein.',
        'split_prodigal_jobs': 'Run multiple prodigal jobs simultaneously for contigs in a fasta file.',
        'include_loose': 'Include loose hits in addition to strict and perfect hits .',
        'include_nudge': 'Include hits nudged from loose to strict hits.',
        'low_quality': 'Use for short contigs to predict partial genes.',
        'num_threads': 'Number of threads (CPUs) to use in the BLAST search.'},
    output_descriptions={'amr_annotations': 'AMR Annotation as .txt and .json file.'},
    name='Annotate MAGs with antimicrobial resistance gene information from CARD.',
    description='Annotate MAGs with antimicrobial resistance gene information from CARD.',
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

#Registrations
plugin.register_semantic_types(CARDDatabase, CARDAnnotation)

plugin.register_semantic_type_to_format(
    CARDDatabase,
    artifact_format=CARDDatabaseDirectoryFormat)
plugin.register_semantic_type_to_format(
    SampleData[CARDAnnotation],
    artifact_format=CARDAnnotationDirectoryFormat)
plugin.register_formats(CARDAnnotationTXTFormat, CARDAnnotationDirectoryFormat,
                        CARDAnnotationJSONFormat,
                        CARDDatabaseFormat, CARDDatabaseDirectoryFormat)

importlib.import_module('q2_amr.types._transformer')
