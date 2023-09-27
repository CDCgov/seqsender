#!/usr/bin/env python3

"""
    Upload consensus sequences and metadata to GISAID's EpiCoV
    Copyright (C) 2022 Freunde von GISAID e.V.
"""

from setuptools import setup, find_packages
import cli3

LONG_DESCRIPTION = open("README.md").read()

setup(

    name=cli3.__name__,
    version=cli3.__version__,
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "cli3 = cli3.__main__:main"
        ]
    },

    description=cli3.__description__,
    long_description=LONG_DESCRIPTION,
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics",
                 "Topic :: Scientific/Engineering :: Medical Science Apps.",
                 "Intended Audience :: Science/Research/Public_Health"],
    keywords=["GISAID",
              "upload",
              "consensus",
              "metadata",
              "hCoV-19",
              "SARS-CoV-2",
              "cli3"],
    author=cli3.__author__,
    author_email=cli3.__author_email__,
    license=cli3.__license__,
    package_data={"": ["*.fa*",
                       "*.csv"]},
    install_requires=["requests==2.28.1",
                      "pandas==1.5.0"],
)

# if __name__ == "__main__":
#     setup()
