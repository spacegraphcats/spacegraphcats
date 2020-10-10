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
    "Programming Language :: Python :: 3.7",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name = 'spacegraphcats',
    version = "1.1",
    description="tools for biological assembly graph neighborhood analysis",
    url="https://github.com/spacegraphcats/spacegraphcats",
    author="C. Titus Brown, Dominik Moritz, Michael P. O'Brien, Felix Reidl, Taylor Reiter, and Blair D. Sullivan",
    author_email="titus@idyll.org,blair.d.sullivan@gmail.com",
    license="BSD 3-clause",
    packages = find_packages(),
    classifiers = CLASSIFIERS,
    entry_points = {'console_scripts': [
        'spacegraphcats  = spacegraphcats.__main__:main'
        ]
    },
    include_package_data=True,
    package_data = { "spacegraphcats": ["Snakefile", "*.yaml", "*.json"] },
    setup_requires = [ "setuptools>=38.6.0",
                       'setuptools_scm', 'setuptools_scm_git_archive' ],
    use_scm_version = {"write_to": "spacegraphcats/version.py"},
    install_requires = [
        'Cython', 'mypy', 'screed', 'pytest',
        'numpy',
        'pandas',
        'snakemake',
        'sortedcontainers',
        'sourmash', 'khmer', 'bbhash']
)
