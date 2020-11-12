import os
from collections import defaultdict
import csv
from typing import Dict, Set, List


class CAtlas:
    """CAtlas class for searching."""

    def __init__(
        self, catlas_directory, load_domfile=True, load_sizefile=False, min_abund=0.0
    ):
        self.name = catlas_directory

        # catlas node ID -> parent
        self.parent = {}  # type: Dict[int, int]

        # catlas node IDs -> { children node IDs }
        self.children = defaultdict(set)  # type: Dict[int, Set[int]]

        # catlas node IDs -> catlas level
        self.levels = {}  # type: Dict[int, int]

        # mapping from cDBG nodes to catlas node IDs, for internal use
        # note: not all cDBG nodes are represented in catlas!
        self._cdbg_to_catlas = {}  # type: Dict[int, int]

        catlas_file = os.path.join(catlas_directory, "catlas.csv")
        self.__load_catlas(catlas_file)
        if load_domfile:
            domfile = os.path.join(catlas_directory, "first_doms.txt")
            self.__load_first_level(domfile)
        if load_sizefile is not None:
            sizefile = os.path.join(catlas_directory, "contigs.info.csv")
            self.__load_size_info(sizefile, min_abund)

    def __load_catlas(self, catlas_file):
        self.max_level = -1
        self.root = -1
        # load everything from the catlas file
        for line in open(catlas_file, "rt"):
            node_id, cdbg_id, level, children = line.strip().split(",")
            # parse out the children
            node_id = int(node_id)
            children = children.strip()
            if children:
                children = children.split(" ")
                children = set(map(int, children))

                # save node -> children, and level
                self.children[node_id] = children
                for child in children:
                    self.parent[child] = node_id
            # since self.children is a defaultdict, we don't have to do
            # anything at the leaves when there are no children

            level = int(level)
            self.levels[node_id] = level
            # update max_node/max_level
            if level > self.max_level:
                self.max_level = level
                self.root = node_id

            # build _cdbg_to_catlas mapping
            if level == 1:
                self._cdbg_to_catlas[int(cdbg_id)] = node_id

    def __load_first_level(self, domfile):
        """
        Load the mapping between first layer catlas and the original DBG nodes.
        """
        # mapping from catlas node IDs to cdbg nodes
        self.layer1_to_cdbg = {}  # type: Dict[int, Set[int]]

        # mapping from cdbg nodes to catlas node IDs
        self.cdbg_to_layer1 = {}  # type: Dict[int, int]

        fp = open(domfile, "rt")
        for line in fp:
            dom_node, *beneath = line.strip().split(" ")

            dom_node = int(dom_node)
            beneath = set(map(int, beneath))

            equiv_cdbg_to_catlas = self._cdbg_to_catlas[dom_node]
            self.layer1_to_cdbg[equiv_cdbg_to_catlas] = beneath
            for cdbg_id in beneath:
                self.cdbg_to_layer1[cdbg_id] = equiv_cdbg_to_catlas

    def __load_size_info(self, sizefile, min_abund):
        kmer_sizes = {}
        weighted_kmer_sizes = {}
        # load size information from file
        with open(sizefile, "rt") as fp:
            reader = csv.DictReader(fp)
            for row in reader:
                contig_id = int(row["contig_id"])
                n_kmers = int(row["n_kmers"])
                mean_abund = float(row["mean_abund"])
                if not min_abund or mean_abund >= min_abund:
                    kmer_sizes[contig_id] = n_kmers
                    weighted_kmer_sizes[contig_id] = mean_abund * n_kmers

        # propagate upwards
        self.kmer_sizes = {}
        self.weighted_kmer_sizes = {}
        for node_id in self:
            level = self.levels[node_id]
            if level == 1:  # aggregate across cDBG nodes
                total_kmers = 0
                total_weighted_kmers = 0
                for cdbg_node in self.layer1_to_cdbg.get(node_id):
                    total_kmers += kmer_sizes.get(cdbg_node, 0)
                    total_weighted_kmers += weighted_kmer_sizes.get(cdbg_node, 0)

                self.kmer_sizes[node_id] = total_kmers
                self.weighted_kmer_sizes[node_id] = total_weighted_kmers
            else:  # aggregate across children
                sub_size = 0
                sub_weighted_size = 0
                for child_id in self.children[node_id]:
                    sub_size += self.kmer_sizes[child_id]
                    sub_weighted_size += self.weighted_kmer_sizes[child_id]
                self.kmer_sizes[node_id] = sub_size
                self.weighted_kmer_sizes[node_id] = sub_weighted_size

    def __iter__(self):
        """
        Iterate through the nodes of the CAtlas such that each node is
        processed after its children.
        """
        Q = [self.root]
        pos = 0
        while pos < len(self.levels):
            v = Q[pos]
            Q.extend(self.children[v])
            pos += 1
        for v in reversed(Q):
            yield v

    def __len__(self):
        return len(self.parent)

    def decorate_with_shadow_sizes(self):
        self.shadow_sizes = {}
        for node_id in self:
            level = self.levels[node_id]
            if level == 1:
                self.shadow_sizes[node_id] = len(self.layer1_to_cdbg[node_id])
            else:
                sub_size = 0
                for child_id in self.children[node_id]:
                    sub_size += self.shadow_sizes[child_id]
                self.shadow_sizes[node_id] = sub_size

    def decorate_with_index_sizes(self, index):
        self.index_sizes = index.build_catlas_node_sizes(self)

    def leaves(self, nodes: List[int] = None) -> Set[int]:
        """
        Return leaves of this CAtlas.
        If nodes is specified, return only those leaves that are descendants of
        the specified nodes; otherwise, return all of them.
        """
        if nodes is None:
            nodes = [self.root]
        leaves = set()  # type: Set[int]
        seen_nodes = set()  # type: Set[int]

        def add_to_shadow(node_id: int):
            if node_id in seen_nodes:
                return
            seen_nodes.add(node_id)

            children_ids = self.children[node_id]

            if len(children_ids) == 0:
                leaves.add(node_id)
            else:
                for child in children_ids:
                    add_to_shadow(child)

        for node in nodes:
            add_to_shadow(node)

        return leaves

    def shadow(self, nodes: List[int]) -> Set[int]:
        """
        Return the cDBG vertices in the shadow of the specified nodes.
        """
        leaves = self.leaves(nodes)
        shadow = set()
        for leaf in leaves:
            shadow.update(self.layer1_to_cdbg[leaf])
        return shadow
