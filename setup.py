from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules=[
    Extension("graph_parser", ["spacegraphcats/graph_parser.pyx"]),
    Extension("walk_dbg", ["spacegraphcats/walk_dbg.pyx"]),
    Extension("graph", ["spacegraphcats/graph.pyx"]),
]

setup(
    name = 'spacegraphcats',
    packages = find_packages(),
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)
