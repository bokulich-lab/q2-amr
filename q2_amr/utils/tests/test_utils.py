from io import StringIO

import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_amr.utils.utils import combine_dataframes, create_df


class MetadataUtilsTest(TestPluginBase):
    package = "q2_amr.utils.tests"

    def setUp(self):
        super().setUp()
        # Setup test data
        self.file_data_1 = "col1\tcol2\nval1\tval2\nval3\tval4"
        self.file_data_2 = "col1\tcol2\nval5\tval6\nval7\tval8"

        self.df1 = pd.read_csv(StringIO(self.file_data_1), sep="\t")
        self.df2 = pd.read_csv(StringIO(self.file_data_2), sep="\t")

        self.df_list = []

        self.df1 = pd.DataFrame(
            {
                "id_name": ["id_value_1", "id_value_1"],
                "col1": ["val1", "val3"],
                "col2": ["val2", "val4"],
            }
        )

        self.df2 = pd.DataFrame(
            {
                "id_name": ["id_value_2", "id_value_2"],
                "col1": ["val5", "val7"],
                "col2": ["val6", "val8"],
            }
        )

        self.expected_combined_df = pd.DataFrame(
            {
                "id_name": ["id_value_1", "id_value_1", "id_value_2", "id_value_2"],
                "col1": ["val1", "val3", "val5", "val7"],
                "col2": ["val2", "val4", "val6", "val8"],
            }
        )

        self.expected_combined_df.index = self.expected_combined_df.index.astype(str)
        self.expected_combined_df.index.name = "id"

    def test_create_df(self):
        # Test create_df function
        create_df(StringIO(self.file_data_1), self.df_list, "id_name", "id_value_1")
        create_df(StringIO(self.file_data_2), self.df_list, "id_name", "id_value_2")

        pd.testing.assert_frame_equal(self.df_list[0], self.df1)
        pd.testing.assert_frame_equal(self.df_list[1], self.df2)

    def test_combine_dataframes(self):
        # Prepare data

        df_list = [self.df1, self.df2]

        combined_df = combine_dataframes(df_list)

        pd.testing.assert_frame_equal(combined_df, self.expected_combined_df)
