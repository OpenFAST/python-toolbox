""""Distribution setup"""


import os

from setuptools import setup

ROOT = os.path.abspath(os.path.dirname(__file__))

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

with open(os.path.join(ROOT, "VERSION")) as version_file:
    VERSION = version_file.read().strip()


setup(
    name="pyfast",
    description="pyfast",
    long_description=LONG_DESCRIPTION,
    version=VERSION,
    url="https://github.com/RHammond2/pyfast/",
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
    packages=["pyfast"],
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "bokeh==1.4.0",
        "pre-commit",
        "black",
        "isort",
        "pytest",
        "pytest-cov",
        "pytest-xdist",
    ],
    test_suite="pytest",
    tests_require=["pytest", "pytest-xdist", "pytest-cov"],
    entry_points={"console_scripts": ["pyfast = pyfast.__main__:main"]},
)
