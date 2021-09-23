# Author: Kobie Kirven
# Davenport Lab - Penn State University
# Date: 9-2-2021

import setuptools
import sys
from pathlib import Path

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


CURRENT_DIR = Path(__file__).parent
sys.path.insert(0, str(CURRENT_DIR))


setuptools.setup(
    name="scims",
    version="0.01",
    packages=["src", "src/scims"],
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "scims=src.__main__:scims",
        ],
    },
    author="Kobie Kirven, Kyle McGovern",
    description="SCiMS: Sex Calling for Metagenomic Sequences",
    install_requires=["setuptools", "biopython", "pandas"],
    python_requires=">=3.5",
)
