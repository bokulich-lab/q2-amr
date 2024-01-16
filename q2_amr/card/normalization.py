from typing import Union

import biom
import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)


def normalize(
    reads: Union[
        SingleLanePerSamplePairedEndFastqDirFmt, SingleLanePerSampleSingleEndFastqDirFmt
    ],
    table: biom.Table,
) -> pd.DataFrame:
    # manifest = reads.manifest.view(pd.DataFrame)
    # #paired = isinstance(reads, SingleLanePerSamplePairedEndFastqDirFmt)
    # for samp in list(manifest.index):
    #     # fwd = manifest.loc[samp, "forward"]
    #     # rev = manifest.loc[samp, "reverse"] if paired else None
    normalized_table = pd.DataFrame(table)
    return normalized_table
