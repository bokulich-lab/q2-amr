# ----------------------------------------------------------------------------
# Copyright (c) 2022, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from q2_types.feature_data import Sequence, FeatureData, ProteinFASTAFormat, DNAFASTAFormat
from qiime2.core.type import Str, Choices, Bool

from q2_amr.card import fetch_card_data, card_annotation
from qiime2.plugin import Citations, Plugin

from q2_amr import __version__

citations = Citations.load("citations.bib", package="q2_amr")

plugin = Plugin(
    name="amr",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-amr",
    package="q2_amr",
    description="This is a QIIME 2 plugin that annotates microbiome sequence data with antimicrobial resistance gene information from CARD.",
    short_description="This is a QIIME 2 plugin that annotates microbiome sequence data with antimicrobial resistance gene information from CARD.",
)

plugin.methods.register_function(
    function=card_annotation,
    inputs={'sequences' : FeatureData[Sequence]},
    parameters={'alignment_tool': Str % Choices(['Blast', 'DIAMOND']),
                'input_type': Str % Choices(['contig', 'protein']),
                'split_prodigal_jobs': Bool,
                'loose': Bool,
                'nudge': Bool,
                'low_quality': Bool,
                'threads': Str},
    outputs=[('rgi_output', pd.DataFrame),
             ('protein_fasta', ProteinFASTAFormat),
             ('dna_fasta', DNAFASTAFormat)],
    input_descriptions={'sequences': 'Sequences to be annotated with rgi.'},
    parameter_descriptions={
                'alignment_tool': 'Specify alignment tool BLAST or DIAMOND (default = BLAST).',
                'input_type': 'Specify data input type contig or protein (default = contig).',
                'split_prodigal_jobs': 'Run multiple prodigal jobs simultaneously for contigs in a fasta file ('
                                       'default: False).',
                'loose': 'Include loose hits in addition to strict and perfect hits (default: False).',
                'nudge': 'Include hits nudged from loose to strict hits (default: False).',
                'low_quality': 'Use for short contigs to predict partial genes (default: False).',
                'threads': 'Number of threads (CPUs) to use in the BLAST search (default=8).'},
    output_descriptions={
        'rgi_output': 'Dataframe with the RGI output.',
        'protein_fasta': 'FASTA file with predicted protein sequences and ORF_ID and ARO accession in the Header.',
        'dna_fasta': 'FASTA file with predicted dna sequences and ORF_ID and ARO accession in the Header.'},
    name='Annotation of sequence data with antimicrobial resistance gene information from CARD.',
    description=('Annotation of sequence data with antimicrobial resistance gene information from CARD.'),
    citations=['Alcock et al. 2023. CARD 2023: expanded curation, support for machine learning, and resistome '
               'prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 51, '
               'D690-D699 [PMID 36263822]']
)