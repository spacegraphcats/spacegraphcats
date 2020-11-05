from setuptools import setup, find_packages

# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: C++",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="spacegraphcats",
    description="tools for biological assembly graph neighborhood analysis",
    url="https://github.com/spacegraphcats/spacegraphcats",
    author="C. Titus Brown, Dominik Moritz, Michael P. O'Brien, Felix Reidl, Taylor Reiter, and Blair D. Sullivan",
    author_email="titus@idyll.org,blair.d.sullivan@gmail.com",
    license="BSD 3-clause",
    packages=find_packages(),
    classifiers=CLASSIFIERS,
    entry_points={
        "console_scripts": ["spacegraphcats  = spacegraphcats.__main__:main"]
    },
    include_package_data=True,
    package_data={"spacegraphcats": ["Snakefile", "*.yaml", "*.json"]},
    setup_requires=[
        "setuptools>=38.6.0",
        "setuptools_scm",
        "setuptools_scm_git_archive",
    ],
    use_scm_version={"write_to": "spacegraphcats/version.py"},
    install_requires=[
        "screed",
        "pytest",
        "numpy",
        "snakemake",
        "sortedcontainers",
        "sourmash",
        "khmer",
        "bbhash >= 0.5",
    ],
    long_description=long_description,
    long_description_content_type="text/markdown",
)
