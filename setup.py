from setuptools import setup, find_packages

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: C++",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name = 'spacegraphcats',
    version = "1.0",
    description="tools for biological assembly graph neighborhood analysis",
    url="https://github.com/spacegraphcats/spacegraphcats",
    author="C. Titus Brown, Dominik Moritz, Michael P. O'Brien, Felix Reidl, Taylor Reiter, and Blair D. Sullivan",
    author_email="titus@idyll.org,blair.d.sullivan@gmail.com",
    license="BSD 3-clause",
    packages = find_packages(),
    classifiers = CLASSIFIERS,
)
