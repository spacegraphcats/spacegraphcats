#! /usr/bin/env python
import sourmash
import pickle
import argparse
from sourmash.search import gather_databases


def main():
    p = argparse.ArgumentParser()
    p.add_argument('node_mh_pickle')
    p.add_argument('sbt')
    args = p.parse_args()

    node_mhs = pickle.load(open(args.node_mh_pickle, 'rb'))
    sbt_obj = sourmash.load_sbt_index(args.sbt)
    databases = ((sbt_obj, args.sbt, 'SBT'),)

    for v in node_mhs.values():
        ss = sourmash.SourmashSignature(v)

        sum_f_uniq = 0.0
        best_match = None
        for (result, _, _, _) in gather_databases(ss, databases, 0, True):
            if best_match is None:
                best_match = result

            sum_f_uniq += result.f_unique_to_query

        print('{} {:.1f}'.format(best_match.name, best_match.f_unique_to_query / sum_f_uniq * 100))


if __name__ == '__main__':
    main()
