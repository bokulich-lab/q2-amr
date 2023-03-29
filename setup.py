# ----------------------------------------------------------------------------
# Copyright (c) 2022, Vinzent Risch.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

setup(
    name="q2-amr",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Vinzent Risch",
    author_email="risch.vinzent@gmail.com",
    description="This is a template for building a new QIIME 2 plugin.",
    url="https://github.com/bokulich-lab/q2-plugin-template",
    entry_points={
        "qiime2.plugins": ["q2-amr=q2_amr.plugin_setup:plugin"]
    },
    package_data={
        "q2_amr": ["citations.bib"],
    },
    zip_safe=False,
)
