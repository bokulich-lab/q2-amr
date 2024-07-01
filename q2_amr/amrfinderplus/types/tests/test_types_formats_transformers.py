# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.plugin.testing import TestPluginBase

from q2_amr.amrfinderplus.types._format import AMRFinderPlusDatabaseDirectoryFormat


class TestAMRFinderPlusDatabaseTypesAndFormats(TestPluginBase):
    package = "q2_amr.amrfinderplus.types.tests"

    def test_amrfinderplus_database_directory_format_validate_positive(self):
        format = AMRFinderPlusDatabaseDirectoryFormat(
            "/Users/rischv/Documents/data/amrfinder/database", mode="r"
        )
        format.validate()
