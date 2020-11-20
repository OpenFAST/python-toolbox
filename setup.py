""""Distribution setup"""


import os

from setuptools import setup

ROOT = os.path.abspath(os.path.dirname(__file__))

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

with open(os.path.join(ROOT, "VERSION")) as version_file:
    VERSION = version_file.read().strip()


setup(
    name="pyFAST",
    description="pyFAST",
    long_description=LONG_DESCRIPTION,
    version=VERSION,
    url="https://github.com/RHammond2/pyFAST/",
    author="Rob Hammond",
    author_email="robert.hammond@nrel.gov",
    classifiers=[
        "Topic :: Utilities",
        "Topic :: Software Development :: Testing",
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Software Development :: Version Control :: Git",
    ],
    packages=["pyFAST"],
    python_requires=">=3.6",
    install_requires=[
        "numpy>=1.15.2",
        "bokeh==1.4.0",
        "future",
        "pandas",
        "matplotlib",
        "chardet",
        "scipy",
        "sympy",
        "pytest"
    ],
    extras_require={
        "dev": ["pre-commit", "black", "isort", "pytest-cov", "pytest-xdist"]
    },
    test_suite="pytest",
    tests_require=["pytest", "pytest-xdist", "pytest-cov"],
    entry_points={"console_scripts": ["pyFAST = pyFAST.__main__:main"]},
)
