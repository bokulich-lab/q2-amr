# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
from pathlib import Path

import qiime2

from q2_amr.amrfinderplus.types import AMRFinderPlusAnnotationsDirFmt
from q2_amr.plugin_setup import plugin
from q2_amr.utils.utils import combine_dataframes, create_append_df


@plugin.register_transformer
def _1(data: AMRFinderPlusAnnotationsDirFmt) -> qiime2.Metadata:
    return qiime2.Metadata(_transfomer_helper(data))


def _transfomer_helper(data):
    df_list = []
    for file_dir_name in os.listdir(str(data)):
        # Check the directory structure
        if os.path.isdir(os.path.join(str(data), file_dir_name)):
            for file in glob.glob(os.path.join(str(data), file_dir_name, "*")):
                file_name = Path(file).stem
                # Annotations file from sample data mags
                if file_name.endswith("_amr_annotations"):
                    id_name = "Sample/MAG name"
                    id_value = file_dir_name + "/" + file_name[:-16]
                # Mutations file from sample data mags
                elif file_name.endswith("_amr_all_mutations"):
                    id_name = "Sample/MAG name"
                    id_value = file_dir_name + "/" + file_name[:-18]
                # Mutations or annotations file from sample data contigs
                else:
                    id_name = "Sample name"
                    id_value = file_dir_name

                create_append_df(
                    file_path=file,
                    df_list=df_list,
                    id_name=id_name,
                    id_value=id_value,
                )
        else:
            # Annotations file from feature data mags
            if file_dir_name.endswith("_amr_annotations.tsv"):
                id_value = file_dir_name[:-20]
            # Mutations file from feature data mags
            else:
                id_value = file_dir_name[:-22]

            create_append_df(
                file_path=os.path.join(str(data), file_dir_name),
                df_list=df_list,
                id_name="MAG name",
                id_value=id_value,
            )

    return combine_dataframes(df_list)
