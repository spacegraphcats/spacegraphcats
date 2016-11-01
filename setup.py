from setuptools import setup, find_packages
from Cython.Build import cythonize

setup(
    name = 'spacegraphcats',
    packages = find_packages(),
    ext_modules = cythonize('cydata.pyx', language='c++')
)
