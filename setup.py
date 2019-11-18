""""Distribution setup"""


import os

from setuptools import setup

ROOT = os.path.abspath(os.path.dirname(__file__))

with open("README.md", "r") as fh:
    long_description = fh.read()

with open(os.path.join(ROOT, "VERSION")) as version_file:
    __version__ = version_file.read().strip()


setup(
    name="pyfast",
    description="pyfast",
    version=__version__,
    url="https://github.com/RHammond2/pyfast/",
    author="Rob Hammond",
    author_email="robert.hammond@nrel.gov",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.7",
    ],
    packages=["pyfast"],
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
