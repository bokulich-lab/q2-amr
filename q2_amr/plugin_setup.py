# ----------------------------------------------------------------------------
# Copyright (c) 2022, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_amr.card import fetch_card_data
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

