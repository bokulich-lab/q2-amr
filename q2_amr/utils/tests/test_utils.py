from io import StringIO
from unittest.mock import call, patch

import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_amr.utils.utils import (
    EXTERNAL_CMD_WARNING,
    colorify,
    combine_dataframes,
    create_count_table,
    create_df,
    read_in_txt,
    run_command,
)


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


class TestAMRAnalysisFeatureTable(TestPluginBase):
    package = "q2_amr.utils.tests"

    @classmethod
    def setUpClass(cls):
        cls.allele_count_df = pd.DataFrame(
            {
                "Reference Sequence": [
                    "ARO:3000796|ID:121|Name:mdtF|NCBI:U00096.1",
                    "ARO:3000815|ID:154|Name:mgrA|NCBI:BA000018.3",
                    "ARO:3000805|ID:172|Name:OprN|NCBI:AE004091.2",
                    "ARO:3000026|ID:377|Name:mepA|NCBI:AY661734.1",
                ],
                "sample1": ["1", "1", "1", "1"],
            }
        )

        cls.gene_count_df = pd.DataFrame(
            {
                "ARO Term": ["mdtF", "mgrA", "OprN", "mepA"],
                "sample1": ["1", "1", "1", "1"],
            }
        )

        cls.mag_count_df = pd.DataFrame(
            {
                "Best_Hit_ARO": ["mdtF", "OprN", "mepA"],
                "sample1/bin1": ["2", "1", "1"],
            }
        )

        cls.frequency_table = pd.DataFrame(
            {
                "sample_id": ["sample1", "sample2"],
                "mdtF": ["1", "0"],
                "mgrA": ["1", "1"],
                "OprN": ["1", "1"],
                "mepA": ["1", "1"],
                "mdtE": ["0", "1"],
            }
        )
        cls.frequency_table.set_index("sample_id", inplace=True)

    def test_read_in_txt_mags(self):
        # Test read_in_txt with output data from annotate_mags_card
        self.read_in_txt_test_body(
            filename="output.mags.txt",
            samp_bin_name="sample1/bin1",
            exp=self.mag_count_df,
            data_type="mags",
            colname="Best_Hit_ARO",
        )

    def test_read_in_txt_reads_allele(self):
        # Test read_in_txt with allele mapping output data from annotate_reads_card
        self.read_in_txt_test_body(
            filename="output.allele_mapping_data.txt",
            samp_bin_name="sample1",
            exp=self.allele_count_df,
            data_type="reads",
            colname="Reference Sequence",
        )

    def test_read_in_txt_reads_gene(self):
        # Test read_in_txt with gene mapping output data from annotate_reads_card
        self.read_in_txt_test_body(
            filename="output.gene_mapping_data.txt",
            samp_bin_name="sample1",
            exp=self.gene_count_df,
            data_type="reads",
            colname="ARO Term",
        )

    def read_in_txt_test_body(self, filename, samp_bin_name, exp, data_type, colname):
        # Create expected and observed count dataframes and compare them
        obs = read_in_txt(
            self.get_data_path(filename), samp_bin_name, data_type, colname
        )
        pd.testing.assert_frame_equal(exp, obs)

    def test_create_count_table(self):
        # Create list of dataframes to be used by create_count_table
        df_list = [self.gene_count_df.copy(), self.gene_count_df.copy()]
        df_list[1].iloc[0, 0] = "mdtE"
        df_list[1].rename(columns={"sample1": "sample2"}, inplace=True)

        # Create observed count table with create_count_table function
        obs = create_count_table(df_list)
        obs = obs.astype(str)

        # Define expected count table
        exp = self.frequency_table

        # Compare expected and observed count table
        pd.testing.assert_frame_equal(exp, obs)

    def test_create_count_table_value_error(self):
        # Assert if ValueError is called when empy list is passed
        self.assertRaises(ValueError, create_count_table, [])


class TestMiscellaneous(TestPluginBase):
    package = "q2_amr.utils.tests"

    def test_colorify(self):
        # Test if colorify function correctly adds color codes
        string = "Hello, world!"
        colored_string = colorify(string)
        expected_output = "\033[1;32mHello, world!\033[0m"
        self.assertEqual(colored_string, expected_output)

    @patch("builtins.print")
    @patch("subprocess.run")
    def test_run_command_verbose(self, mock_subprocess_run, mock_print):
        cmd = ["echo", "hello"]

        run_command(cmd, verbose=True)

        # Check that the warning and the command were printed
        expected_calls = [
            call(EXTERNAL_CMD_WARNING),
            call("\nCommand:", end=" "),
            call("echo hello", end="\n\n"),
        ]
        mock_print.assert_has_calls(expected_calls)

        # Check that subprocess.run was called with the correct parameters
        mock_subprocess_run.assert_called_once_with(cmd, check=True, cwd=None)

    @patch("subprocess.run")
    def test_run_command_with_cwd(self, mock_subprocess_run):
        cmd = ["echo", "hello"]
        cwd = "/some/directory"

        run_command(cmd, cwd=cwd, verbose=False)

        # Check that subprocess.run was called with the correct parameters including cwd
        mock_subprocess_run.assert_called_once_with(cmd, check=True, cwd=cwd)
