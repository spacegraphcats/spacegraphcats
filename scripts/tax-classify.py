#! /usr/bin/env python
"""
Calculate the "taxonomic purity" of bins built by characterize_catlas_regions.
"""
import screed
import argparse
import sourmash_lib
import numpy
from itertools import product
from sourmash_lib._minhash import hash_murmur
from sourmash_lib.lca import lca_utils
from pickle import dump
import pickle
import collections

taxonomy_levels = ['superkingdom', 'phylum', 'class', 'order', 'family',
                   'genus', 'species', 'strain']


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('lca_db')
    args = parser.parse_args()

    # load the minhashes calculated by search.characterize_catlas_regions
    group_ident = pickle.load(open(args.filename + '.node_mh', 'rb'))

    # load the LCA database
    dblist, ksize, scaled = lca_utils.load_databases([args.lca_db], None)
    db = dblist[0]

    n_empty = 0
    pure_at_rank = collections.defaultdict(int)

    # for each group,
    for g in group_ident.values():
        # collect all the hash values for this partition
        mins = g.get_mins()

        # now, collate the lineage assignments across all the hashes
        assignments = collections.Counter()
        for m in mins:
            x = db.get_lineage_assignments(m)
            for lineage in x:
                assignments[lineage] += 1

        if not assignments:
            # empty, no hash assignments for this group!
            n_empty += 1
            continue

        # for each taxonomic level, calculate purity
        for level, rank in enumerate(reversed(taxonomy_levels)):
            index = len(taxonomy_levels) - level

            # track number of different items at this tax rank
            level_counter = collections.Counter()

            count = 0
            for lineage, count in assignments.items():
                # trim the lineage back to the given rank
                sublineage = lineage[:index]
                assert sublineage[-1].rank == rank, lineage

                # count the trimmed lineage!
                level_counter[sublineage] += count

            ### Note for @sgsmob: at this point, 'level_counter' contains
            ### all lineages trimmed back to 'rank', with associated counts.
            ### This is a collections.Counter object, so both .items() and
            ### the method '.most_common()' will return lists of
            ### (lineage, count), with 'most_common()' returning them in
            ### highest-to-lowest order.

            ### Here is some code to compute the fraction of the bin that
            ### belongs to the top-ranked lineage:

            top_lineage, top_count = next(iter(level_counter.most_common()))
            total = sum(level_counter.values())
            lineage_display = "; ".join(lca_utils.zip_lineage(top_lineage))
            #print(rank, '{:.1f}%'.format(top_count / total * 100), lineage_display)

            if len(level_counter) == 0:   # should never get here!
                assert 0, assignments

            # if we have a 100% pure bin, quit now.
            if len(level_counter) == 1:
                pure_at_rank[rank] += 1

                # print out the egregiously bad ones just to double check...
                if rank in ('phylum', 'superkingdom'):
                    print('---')
                    for assignment in assignments:
                        print('\t', '; '.join(lca_utils.zip_lineage(assignment)))
                break

    # calculate summary numbers!
    total = len(group_ident)              # total number of bins
    print('total bins: {}'.format(total))
    print('unassigned bins: {}'.format(n_empty))

    sum_pure = 0
    for rank in reversed(taxonomy_levels):
        print('pure at rank {}: {}'.format(rank, pure_at_rank[rank]))
        sum_pure += pure_at_rank[rank]

    print('not pure at any rank:', total - sum_pure - n_empty)


if __name__ == '__main__':
    main()
