from __future__ import print_function

import sys
import khmer, khmer.utils
import screed
from collections import OrderedDict, defaultdict
import os, os.path
from spacegraphcats.catlas.graph_parser import write
import gzip


class Pathfinder(object):
    "Track segment IDs, adjacency lists, and assembled contigs."
    def __init__(self, ksize, gxtfile, contigfile=None, assemble=False):
        self.ksize = ksize

        self.node_counter = 0
        self.nodes = {}                      # node IDs (int) to size
        self.nodes_to_kmers = {}             # node IDs (int) to kmers
        self.kmers_to_nodes = {}             # kmers to node IDs
        self.adjacencies = defaultdict(set)  # node to node
        self.labels = defaultdict(set)       # nodes to set of labels
        self.assemblyfp = None
        if assemble:
            self.assemblyfp = gzip.open(contigfile, 'wt')
        self.adjfp = open(gxtfile + '.adj', 'wt')

    def new_hdn(self, kmer):
        "Add a new high-degree node to the cDBG."
        if kmer in self.kmers_to_nodes:
            return self.kmers_to_nodes[kmer]

        this_id = self.node_counter
        self.node_counter += 1

        self.nodes[this_id] = 1
        self.nodes_to_kmers[this_id] = kmer
        self.kmers_to_nodes[kmer] = this_id



        return this_id

    def new_linear_node(self):
        "Add a new linear path to the cDBG."
        node_id = self.node_counter
        self.node_counter += 1

        return node_id

    def add_adjacency(self, node_id, adj):
        "Add an edge between two nodes to the cDBG."
        node_id, adj = min(node_id, adj), max(node_id, adj)

        self.adjfp.write('{},{}\n'.format(node_id, adj))

    def add_label(self, kmer, label):
        x = self.labels[kmer]
        x.add(label)

    def add_assembly(self, path_id, contig):
        self.assemblyfp.write('>{}\n{}\n'.format(path_id, contig))


def traverse_and_mark_linear_paths(graph, nk, stop_bf, pathy, degree_nodes):
    size, adj_kmers, visited = graph.traverse_linear_path(nk, degree_nodes,
                                                          stop_bf)
    if not size:                          # 0 length paths
        return

    # get an ID for the new path
    path_id = pathy.new_linear_node()

    # output a contig if requested
    if pathy.assemblyfp:
        asm = khmer.LinearAssembler(graph, stop_bf)
        contig = asm.assemble(nk)
        pathy.add_assembly(path_id, contig)
        if size and not contig:
            print('nonzero size, but contig is not produced. WTF.')
        if contig and graph.get_min_count(contig) == 0:
            print('generated k-mers not in BF. sigh.')
        if contig:
            stop_bf.add(contig[:graph.ksize()])
            stop_bf.add(contig[-graph.ksize():])

        if len(contig) and size + graph.ksize() - 1 != len(contig):
            print('visited k-mers != contig size. WTF?')
    else:
        assert 0

    # add all adjacencies, if any
    if adj_kmers:
        for kmer in adj_kmers:
            adj_node_id = pathy.kmers_to_nodes[kmer]
            pathy.add_adjacency(path_id, adj_node_id)
    # a purely linear path; add all the k-mers to stop bf to prevent
    # traversing in both directions.
    else:
        for k in visited:
            stop_bf.add(k)
        # note: if you knew which k-mer was the other end,
        # you could just add that; but we don't know.



def run(args):

    # @CTB this is kind of a hack - nothing tricky going on, just want to
    # specify memory on the command line rather than graph size...
    graph_tablesize = int(args.memory * 8.0 / 4.0)

    assert args.ksize % 2, "ksize must be odd"

    if args.label:
        label_list = []

    output_dir = args.output
    if not output_dir:
        if len(args.seqfiles) > 1:
            print('** please specify an output directory with -o',
                  file=sys.stderr)
            sys.exit(-1)

        output_dir = os.path.basename(args.seqfiles[0])
        if output_dir.endswith('.fa'):
            output_dir = output_dir[:-3]
        elif output_dir.endswith('.fa.gz'):
            output_dir = output_dir[:-6]

    # set this so we can read it for logging
    args.output = output_dir
    # gxtfile = os.path.basename(output_dir) + '.gxt'
    gxtfile = os.path.join(output_dir, "cdbg.gxt")
    contigfile = os.path.join(output_dir, "contigs.fa.gz")

    print('')
    print('placing output in directory:', output_dir)
    print('gxt will be:', gxtfile)
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        print('(note: directory already exists)')

    print('')
    if args.loadgraph:
        print('loading nodegraph from:', args.loadgraph)
        graph = khmer.Nodegraph.load(args.loadgraph)
        print('creating accompanying stopgraph')
        ksize = graph.ksize()
        hashsizes = graph.hashsizes()
        stop_bf = khmer.Nodegraph(ksize, 1, 1, primes=hashsizes)
    else:
        print('building graphs and loading files')

        # Create graph and a stop bloom filter - one for loading, one for
        # traversing. Create them all here so that we can error out quickly
        # if memory is a problem.

        # @CTB note that hardcoding '2' here is not nec a great idea.
        graph = khmer.Nodegraph(args.ksize, graph_tablesize, 2)
        stop_bf = khmer.Nodegraph(args.ksize, graph_tablesize, 2)
        n = 0

        # load in all of the input sequences, one file at a time.
        for seqfile in args.seqfiles:
            fp = screed.open(seqfile)
            for record in khmer.utils.clean_input_reads(fp):
                if len(record.cleaned_seq) < graph.ksize(): continue
                n += 1
                if n % 100000 == 0:
                    print('...', seqfile, n)
                graph.consume(record.cleaned_seq)
            fp.close()

        # complain if too small set of graphs was used.
        fp_rate = khmer.calc_expected_collisions(graph,
                                                 args.force, max_false_pos=.05)

    ksize = graph.ksize()

    # initialize the object that will track information for us.
    pathy = Pathfinder(ksize, gxtfile, contigfile, not args.no_assemble)

    print('finding high degree nodes')
    if args.label:
        print('(and labeling them, per request)')
    degree_nodes = khmer.HashSet(ksize)
    linear_starts = khmer.HashSet(ksize)
    n = 0
    skipped = 0
    for seqfile in args.seqfiles:
        fp = screed.open(seqfile)
        for record in khmer.utils.clean_input_reads(fp):
            if len(record.cleaned_seq) < ksize:
                skipped += 1
                continue
            n += 1
            if n % 100000 == 0:
                print('...2', seqfile, n)
            # walk across sequences, find all high degree nodes,
            # name them and cherish them.
            these_hdn = graph.find_high_degree_nodes(record.cleaned_seq)
            if these_hdn:
                degree_nodes += these_hdn
            else:
                # possible linear node? check first and last k-mer.
                # (the logic here is that every purely linear node must
                # start or end in *some* record.sequence - so where we have
                # record sequences that have only 1 neighbor, those will be
                # all possible linear nodes.
                first_kmer = record.sequence[:ksize]
                last_kmer = record.sequence[-ksize:]
                assert len(last_kmer) == ksize

                if len(graph.neighbors(first_kmer)) == 1:
                    linear_starts.add(graph.hash(first_kmer))
                if len(graph.neighbors(last_kmer)) == 1:
                    linear_starts.add(graph.hash(last_kmer))

            if args.label:
                label_list.append(record.name)
                for kmer in these_hdn:
                    pathy.add_label(kmer, n)
        fp.close()

    print('read {}, skipped {} for being too short'.format(n, skipped))

    # get all of the degree > 2 kmers and give them IDs.
    for kmer in degree_nodes:
        pathy.new_hdn(kmer)
        stop_bf.add(kmer)

    print('traversing linear segments from', len(degree_nodes), 'nodes')

    # now traverse from each high degree node into all neighboring nodes,
    # seeking adjacencies.  if neighbor is high degree node, add it to
    # adjacencies; if neighbor is not, then traverse the linear path &
    # assemble if desired.
    for n, k in enumerate(degree_nodes):
        if n % 10000 == 0:
            print('...', n, 'of', len(degree_nodes))

        # retrieve the node ID of the primary segment.
        k_id = pathy.kmers_to_nodes[k]

        # here is where we would output this k-mer to the contig file if we
        # wanted to.
        nk_id = pathy.kmers_to_nodes[k]
        k_str = khmer.reverse_hash(k, ksize)
        pathy.add_assembly(nk_id, k_str)

        # find all the neighbors of this high-degree node.
        nbh = graph.neighbors(k)
        for nk in nbh:
            # neighbor is high degree? fine, mark its adjacencies.
            if nk in degree_nodes:
                nk_id = pathy.kmers_to_nodes[nk.kmer_u]
                pathy.add_adjacency(k_id, nk_id)
            else:
                # linear! walk it.
                traverse_and_mark_linear_paths(graph, nk, stop_bf, pathy,
                                               degree_nodes)

    # now, clean up at the end -- make sure we've hit all the possible
    # linear nodes.
    print('traversing from {} potential linear starts'.format(len(linear_starts)))
    for n, k in enumerate(linear_starts):
        traverse_and_mark_linear_paths(graph, k, stop_bf, pathy, degree_nodes)

    print('{} linear segments and {} high-degree nodes'.\
              format(pathy.node_counter, len(pathy.nodes)))

    del graph
    del stop_bf

    # save to GXT.
    print('saving gxtfile', gxtfile)

    all_labels = set()
    label_counts = {}

    pathy.adjfp.close()
    adj_fp = open(gxtfile + '.adj', 'rt')

    # this uniqifies the edges...
    for line in adj_fp:
        a, b = line.split(',')
        a = int(a)
        b = int(b)
        pathy.adjacencies[a].add(b)

    adj_fp.close()
    try:
        os.unlink(gxtfile + '.adj')
    except:
        print('cannot remove', gxtfile + '.adj')

    # ...and now print them out.
    edges = []
    for k, v in pathy.adjacencies.items():
        for dest in v:
            # don't add loops
            if (k != dest):
                edges.append((k, dest))

    with open(gxtfile, 'wt') as fp:
        write(fp, pathy.node_counter, edges)

    if not args.no_assemble:
        pathy.assemblyfp.close()

    if args.label:
        print('note: used/assigned %d labels total' % (len(set(all_labels)),))
        print('counts:', label_counts)

        assert label_list
        print('dumping label list now.')
        label_file = os.path.basename(output_dir) + '.labels.txt'
        label_file = os.path.join(output_dir, label_file)

        with open(label_file, "wt") as fp:
            for n, label in enumerate(label_list):
                fp.write("{} {}\n".format(n + 0, label))
