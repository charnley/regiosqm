#!/usr/bin/env python

import setuptools

setuptools.setup(
    name="regiosqm",
    maintainer="Jimmy Kromann",
    python_requires=">=3.9",
    packages=setuptools.find_packages("src"),
    package_dir={"": "src"},
)
