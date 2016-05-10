#! /usr/bin/env python3
import argparse
import sys, os
from os import path
from graph import Graph
from parser import parse_minhash

DEBUG = True

class AttributeDict(dict): 
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__

def report(msg):
	print(msg)

def debug(msg):
	if DEBUG == True:
		print(msg)

def warn(msg):
	print("  Warning:",msg)

def error(msg):
	raise RuntimeError(msg)

def parse_minhash_dict(file):
	res = {}
	def add_minhash(v, hashlist):
		res[int(v)] = hashlist
	parse_minhash(file, add_minhash)
	return res

def read_project_file(projectpath, filename):
	fullpath = path.join(projectpath, filename)
	if not path.exists(fullpath):
		error("Missing file {} in {}".format(filename, projectpath))
	report("Found {} in {}".format(filename, projectpath))
	return open(fullpath, 'r')

parser = argparse.ArgumentParser()
parser.add_argument('project', help='Project directory. Must contain a .gxt and an .mxt file with the same base name \
									 as the directory.', type=str)
parser.add_argument('n', help="The target atlas size.", type=int )
args = parser.parse_args()

project = AttributeDict()

if not path.isdir(args.project):
	error("{} is not a valid pathname".format(args.project))
if not path.exists(args.project):
	error("Project directory {} does not exist".format(args.project))

project.path = args.project
project.name = path.basename(path.normpath(args.project))

report("Project {} in {}".format(project.name, project.path))

""" Make sure .gxt and .mxt with the right naming conventions exist and load them """

file = read_project_file(project.path, project.name+".gxt")
project.graph, project.node_attr, project.edge_attr = Graph.from_gxt(file)

report("Loaded graph with {} vertices and {} edges".format(len(project.graph),project.graph.num_edges()))

file = read_project_file(project.path, project.name+".mxt")
project.minhashes = parse_minhash_dict(file)

for v in project.graph:
	if v not in project.minhashes:
		warn("Vertex {} is missing minhashes".format(v))

report("Loaded minhashes for graph")

