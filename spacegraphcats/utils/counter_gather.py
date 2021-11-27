from collections import Counter

class CounterGather:
    """
    A refactoring of sourmash.index.CounterGather to use sets of hashes.

    The public interface is `peek(...)` and `consume(...)` only.
    """
    def __init__(self, query_set):
        # track query
        self.orig_query_set = set(query_set)

        # track matches in a list & their locations
        self.setlist = []
        self.locations = []

        # ...and overlaps with query
        self.counter = Counter()

        # cannot add matches once query has started.
        self.query_started = 0

    def add(self, hashes, ident, require_overlap=True):
        "Add this set of hashes in as a potential match."
        if self.query_started:
            raise ValueError("cannot add more hashes to counter after peek/consume")

        # upon insertion, count & track overlap with the specific query.
        overlap = self.orig_query_set & hashes
        if overlap:
            i = len(self.setlist)

            self.counter[i] = len(overlap)
            self.setlist.append(hashes)
            self.locations.append(ident)
        elif require_overlap:
            raise ValueError("no overlap between query and signature!?")

    def peek(self, cur_query_set):
        "Get next 'gather' result for this database, w/o changing counters."
        self.query_started = 1

        # empty? nothing to search.
        counter = self.counter
        if not counter:
            return []

        setlist = self.setlist
        assert setlist

        if not cur_query_set:             # empty query? quit.
            return []

        if cur_query_set & self.orig_query_set != cur_query_set:
            raise ValueError("current query not a subset of original query")

        # Find the best match -
        most_common = counter.most_common()
        dataset_id, match_size = most_common[0]

        # pull match and location.
        match_set = setlist[dataset_id]

        # calculate containment
        cont = len(cur_query_set & match_set) / len(match_set)
        assert cont

        # calculate intersection of this "best match" with query.
        intersect_set = cur_query_set & match_set
        location = self.locations[dataset_id]

        # build result & return intersection
        return cont, match_set, location, intersect_set

    def consume(self, intersect_set):
        "Remove the given hashes from this counter."
        self.query_started = 1

        if not intersect_set:
            return

        setlist = self.setlist
        counter = self.counter

        most_common = counter.most_common()

        # Prepare counter for finding the next match by decrementing
        # all hashes found in the current match in other datasets;
        # remove empty datasets from counter, too.
        for (dataset_id, _) in most_common:
            remaining_set = setlist[dataset_id]
            intersect_count = len(intersect_set & remaining_set)
            if intersect_count:
                counter[dataset_id] -= intersect_count
                if counter[dataset_id] == 0:
                    del counter[dataset_id]

    def do_full_gather(self):
        "Get all matching names for the original query."
        matching_kmers = set(self.orig_query_set)
        filtered_names = {}

        x = self.peek(matching_kmers)
        while x:
            cont, match_set, name, intersect_set = x
            self.consume(intersect_set)
            matching_kmers -= intersect_set
            assert name not in filtered_names
            filtered_names[name] = len(intersect_set)

            x = self.peek(matching_kmers)

        return filtered_names
