#! /usr/bin/env python3
"""
"""
from __future__ import print_function
import argparse
from spacegraphcats import graph_parser
from spacegraphcats.catlas_reader import CAtlasReader
from spacegraphcats.catlas import CAtlas
from spacegraphcats.graph import VertexDict
from collections import defaultdict
import os
import sys

p = argparse.ArgumentParser()
p.add_argument('catlas_prefix', help='catlas prefix')
p.add_argument('catlas_r', type=int, help='catlas radius to load')
args = p.parse_args()   

_radius = args.catlas_r
_basename = os.path.basename(args.catlas_prefix)
_catgxt = '%s.catlas.%d.gxt' % (_basename, _radius)
_catmxt = '%s.catlas.%d.mxt' % (_basename, _radius)
_catgxt = os.path.join(args.catlas_prefix, _catgxt)
_catmxt = os.path.join(args.catlas_prefix, _catmxt)
_catlas = CAtlas.read(_catgxt, _catmxt, _radius)    

s = "{} (@{}): size={}, vertex={}, mh={}"
for n in _catlas.nodes():
    print(s.format(n.id, n.level, n.size, n.vertex, n.minhash))
