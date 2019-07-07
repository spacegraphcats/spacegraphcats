#! /usr/bin/env python
import sourmash
import pickle
import argparse
from sourmash.search import gather_databases
from sourmash.lca.lca_utils import LCA_Database


def main():
    p = argparse.ArgumentParser()
    p.add_argument('node_mh_pickle')
    p.add_argument('lca_db')
    args = p.parse_args()

    node_mhs = pickle.load(open(args.node_mh_pickle, 'rb'))
    lca_obj = LCA_Database()
    lca_obj.load(args.lca_db)
    databases = ((lca_obj, args.lca_db, 'LCA'),)

    for k, v in node_mhs.items():
        ss = sourmash.SourmashSignature(v)

        sum_f_uniq = 0.0
        best_match = None
        for (result, _, _, _) in gather_databases(ss, databases, 0, True):
            if best_match is None:
                best_match = result

            sum_f_uniq += result.f_unique_to_query

        if not best_match:
            print('** no match for {}'.format(k))
            continue

        print('node {} matches {} {:.1f}'.format(k, best_match.name[:30], best_match.f_unique_to_query / sum_f_uniq * 100))


if __name__ == '__main__':
    main()
