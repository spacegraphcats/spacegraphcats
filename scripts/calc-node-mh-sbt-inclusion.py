#! /usr/bin/env python
import sourmash
import pickle
import argparse
from sourmash.search import gather_databases
from sourmash.lca import lca_utils
from sourmash.lca.lca_utils import LCA_Database


def main():
    p = argparse.ArgumentParser()
    p.add_argument("node_mh_pickle")
    p.add_argument("lca_db")
    args = p.parse_args()

    node_mhs = pickle.load(open(args.node_mh_pickle, "rb"))
    lca_obj = LCA_Database()
    lca_obj.load(args.lca_db)
    databases = ((lca_obj, args.lca_db, "LCA"),)

    d = {}
    n_pure95 = 0
    total = 0

    for k, v in node_mhs.items():
        ss = sourmash.SourmashSignature(v)

        results = [x[0] for x in gather_databases(ss, databases, 0, True)]
        sum_f_uniq = sum([result.f_unique_to_query for result in results])

        keep_results = []
        for result in results:
            if result.f_unique_to_query < 0.10:
                break
            keep_results.append(result)

        if not keep_results:
            print("** no match for {}".format(k))
            continue

        idents = [result.name.split()[0].split(".")[0] for result in keep_results]
        idxlist = [lca_obj.ident_to_idx[ident] for ident in idents]
        lidlist = [lca_obj.idx_to_lid[idx] for idx in idxlist]
        lineages = [lca_obj.lid_to_lineage[lid] for lid in lidlist]

        tree = lca_utils.build_tree(lineages)
        lca, reason = lca_utils.find_lca(tree)

        level = "*none*"
        if lca:
            level = lca[-1].rank

        lineage = ";".join(lca_utils.zip_lineage(lca, truncate_empty=True))

        this_f_uniq = sum([result.f_unique_to_query for result in keep_results])

        print(
            "node {} matches {} @ {:.1f}".format(
                k, level, this_f_uniq / sum_f_uniq * 100
            )
        )

        if level in ("strain", "genus", "species") and this_f_uniq / sum_f_uniq >= 0.95:
            n_pure95 += 1
        total += 1

    print("XXX", n_pure95, total)


if __name__ == "__main__":
    main()
