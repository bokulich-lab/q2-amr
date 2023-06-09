# q2-amr
![CI](https://github.com/bokulich-lab/q2-assembly/actions/workflows/ci.yml/badge.svg)
[![codecov](https://codecov.io/gh/bokulich-lab/q2-assembly/branch/main/graph/badge.svg?token=THMBOFUZR0)](https://codecov.io/gh/bokulich-lab/q2-assembly)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

QIIME 2 plugin for antimicrobial resistance gene annotation of MAGs and metagenomic reads.

## Installation
To install _q2-amr_, follow the installation steps described below.

```shell
mamba create -yn q2-amr \
   -c conda-forge -c bioconda -c qiime2 -c defaults -c https://packages.qiime2.org/qiime2/2023.5/tested --strict-channel-priority qiime2 q2cli q2templates q2-types q2-types-genomics q2-amr rgi altair click==7.1.2

mamba activate q2-amr

pip install git+https://github.com/misialq/rgi.git@py38-fix
```

Refresh cache and check that everything worked:
```shell
qiime dev refresh-cache
qiime info
```

## Functionality
This QIIME 2 plugin contains actions used to annotate short single/paired-end
sequencing reads and MAGs with antimicrobial resistance genes. Currently, the [CARD](https://card.mcmaster.ca) database is supported  (for details on
the implementation and usage, please refer to the [rgi](https://github.com/arpcard/rgi) documentation). Below you will
find an overview of actions available in the plugin.

| Action                     | Description                                                                          | Underlying tool                       |
|----------------------------|--------------------------------------------------------------------------------------|---------------------------------------|
| fetch-card-db              | Download CARD data.                                                                  | [rgi](https://github.com/arpcard/rgi) |
| annotate-mags-card         | Annotate MAGs with antimicrobial resistance gene information from CARD.              | [rgi](https://github.com/arpcard/rgi) |
| annotate-reads-card        | Annotate metagenomic reads with antimicrobial resistance gene information from CARD. | [rgi](https://github.com/arpcard/rgi) |
| heatmap                    | Create a heatmap from annotate-mags-card output files.                               |                                       |
| visualize-annotation-stats | Plot annotate-reads-card annotation statistics.                                      |                                       |


## Dev environment
This repository follows the _black_ code style. To make the development slightly easier
there are a couple of pre-commit hooks included here that will ensure that your changes
follow that formatting style. Before you start working on the code, please
install the hooks by executing `make dev` in your conda environment. From then on,
they will be run automatically every time you commit any changes.
