#! /usr/bin/env python3

import argparse, sys, os, re
import gzip, glob
from operator import itemgetter
from os import path

from spacegraphcats.graph import Graph
from spacegraphcats.graph_io import read_from_gxt, write_to_gxt
from spacegraphcats.rdomset import low_degree_orientation, alt_ldo
from spacegraphcats.graph_parser import parse

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

def read_project_file(projectpath, filename):
    """
        Attempts to read a project file. Checks whether compressed file (.gz)
        exists that otherwise matches the filename.
    """
    fullpath = path.join(projectpath, filename)
    zipped = False
    if not path.exists(fullpath):
        if not path.exists(fullpath+".gz"):
            error("Missing file {} in {}".format(filename, projectpath))
        else:
            zipped = True
            fullpath += ".gz"
            filename += ".gz" # for consistent report
    report("Found {} in {}".format(filename, projectpath))
    if zipped:
        return gzip.open(fullpath, 'rt')    
    return open(fullpath, 'r')

def load_and_compute_augg(project):
    """ 
        Returns a project.radius-dtf augmentation of project.graph.
        Loads cached augmentations from the project directory and writes
        newly computed augmentations into it.
    """
    augname = path.join(project.path, project.name + ".aug.{}.ext")

    augs = {}
    for f in glob.glob(augname.format("[0-9]*")):
        d = int(f.split(".")[-2])
        augs[d] = f

    aug0 = augname.format('0')
    if 0 in augs:
        report('  Load cached augmentation from {}'.format(aug0))
        with open(aug0, 'r') as f:
            parse(f, None, project.graph.add_arc)
    else:
        report('  Running low degree orientation and cache it as {}'.format(aug0))
        alt_ldo(project.graph)
        with open(aug0, 'w') as f:
            write_to_gxt(f, project.graph, 1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('project', help='Project directory. Must contain a .gxt and an .mxt file with the same base name \
                                         as the directory.', type=str)
    parser.add_argument('r', help="The catlas' radius.", type=int )
    args = parser.parse_args()

    project = AttributeDict()
    project.radius = args.r

    if not path.isdir(args.project):
        error("{} is not a valid pathname".format(args.project))
    if not path.exists(args.project):
        error("Project directory {} does not exist".format(args.project))

    project.path = args.project
    project.name = path.basename(path.normpath(args.project))

    report("Project {} in {}".format(project.name, project.path))

    """ 
        Make sure .gxt with the right naming conventions exist and load them 
    """

    file = read_project_file(project.path, project.name+".gxt")
    project.graph = read_from_gxt(file, project.radius, False)

    report("Loaded graph with {} vertices, {} arcs and {} components".format(
        len(project.graph),
        project.graph.num_arcs(),
        project.graph.num_components()))


    """ Compute / load r-dominating set """

    report("\nDomset computation\n")

    load_and_compute_augg(project)

    """ Compute catlas """

    report("\nCatlas computation\n")

    sys.exit(0)


if __name__ == '__main__':
    main()
