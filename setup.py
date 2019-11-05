""""Distribution setup"""

from setuptools import setup

import pyFAST

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name="pyFAST",
    version=pyFAST.__version__,
    description="pyFAST",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.7",
    ],
    packages=["pyFAST"],
    install_requires=[
        "numpy",
        "bokeh",
        "pre-commit",
        "black",
        "isort",
        "pytest",
        "pytest-cov",
        "pytest-xdist",
    ],
    test_suite="pytest",
    tests_require=["pytest", "pytest-xdist", "pytest-cov"],
)
