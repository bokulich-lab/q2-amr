# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_types.sample_data import SampleData
from qiime2.plugin import SemanticType

CARDDatabase = SemanticType("CARDDatabase")
CARDKmerDatabase = SemanticType("CARDKmerDatabase")
CARDAnnotation = SemanticType("CARDAnnotation", variant_of=SampleData.field["type"])
CARDAlleleAnnotation = SemanticType(
    "CARDAlleleAnnotation", variant_of=SampleData.field["type"]
)
CARDGeneAnnotation = SemanticType(
    "CARDGeneAnnotation", variant_of=SampleData.field["type"]
)
