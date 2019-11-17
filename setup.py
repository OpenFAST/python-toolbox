""""Distribution setup"""

from setuptools import setup

# import pyfast

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name="pyfast",
    version="0.1.0",
    description="pyfast",
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
