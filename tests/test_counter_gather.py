# tests adapted from sourmash test_index.py
import pytest

from spacegraphcats.utils.counter_gather import CounterGather


def _consume_all(query, counter):
    results = []
    last_intersect_size = None
    while 1:
        result = counter.peek(query)
        if not result:
            break

        cont, match_set, location, intersect_set = result
        print(location, len(intersect_set))
        if last_intersect_size:
            assert len(intersect_set) <= last_intersect_size

        last_intersect_size = len(intersect_set)

        counter.consume(intersect_set)
        query -= intersect_set

        results.append((cont, match_set, location, len(intersect_set)))

    return results


def test_counter_gather_1():
    # check a contrived set of non-overlapping gather results,
    # generated via CounterGather
    query = set(range(0, 20))
    set1 = set(range(0, 10))
    set2 = set(range(10, 15))
    set3 = set(range(15, 17))

    counter = CounterGather(query)
    counter.add(set1, 'match1')
    counter.add(set2, 'match2')
    counter.add(set3, 'match3')

    results = _consume_all(query, counter)
    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results


def test_counter_gather_1_b():
    # check a contrived set of somewhat-overlapping gather results,
    # generated via CounterGather. Here the overlaps are structured
    # so that the gather results are the same as those in
    # test_counter_gather_1(), even though the overlaps themselves are
    # larger.

    query = set(range(0, 20))
    set1 = set(range(0, 10))
    set2 = set(range(7, 15))
    set3 = set(range(13, 17))
    
    # load up the counter
    counter = CounterGather(query)
    counter.add(set1, 'match1')
    counter.add(set2, 'match2')
    counter.add(set3, 'match3')

    results = _consume_all(query, counter)

    expected = (['match1', 10],
                ['match2', 5],
                ['match3', 2],)
    assert len(results) == len(expected), results


def test_counter_gather_exact_match():
    # query == match
    query = set(range(0, 20))
    
    # load up the counter
    counter = CounterGather(query)
    counter.add(query, "somewhere over the rainbow")

    results = _consume_all(query, counter)
    assert len(results) == 1
    (cont, match_set, ident, intersect_mh) = results[0]

    assert ident == 'somewhere over the rainbow'


def test_counter_gather_add_after_peek():
    # query == match
    query = set(range(0, 20))
    
    # load up the counter
    counter = CounterGather(query)
    counter.add(query, "somewhere over the rainbow")

    # cannot add after peek or consume
    counter.peek(query)

    with pytest.raises(ValueError):
        counter.add(query, "try again")


def test_counter_gather_add_after_consume():
    # cannot add after peek or consume
    query = set(range(0, 20))

    # load up the counter
    counter = CounterGather(query)
    counter.add(query, 'somewhere over the rainbow')

    counter.consume(query)

    with pytest.raises(ValueError):
        counter.add(query, "try again")


def test_counter_gather_consume_empty_intersect():
    query = set(range(0, 20))

    # load up the counter
    counter = CounterGather(query)
    counter.add(query, 'somewhere over the rainbow')

    # nothing really happens here :laugh:, just making sure there's no error
    counter.consume(set())


def test_counter_gather_empty_initial_query():
    # check empty initial query
    query = set()

    match = set(range(0, 10))

    # load up the counter
    counter = CounterGather(query)
    counter.add(match, "match1", require_overlap=False)

    assert counter.peek(query) == []


def test_counter_gather_empty_cur_query():
    # test empty cur query
    query = set(range(0, 20))

    # load up the counter
    counter = CounterGather(query)
    counter.add(query, 'somewhere over the rainbow')

    cur_query = set()
    results = _consume_all(cur_query, counter)
    assert results == []


def test_counter_gather_bad_cur_query():
    # test cur query that is not subset of original query
    query = set(range(0, 20))

    # load up the counter
    counter = CounterGather(query)
    counter.add(query, 'somewhere over the rainbow')

    cur_query_set = set(range(20, 30))
    with pytest.raises(ValueError):
        counter.peek(cur_query_set)


def test_counter_gather_add_no_overlap():
    # check adding match with no overlap w/query
    query = set(range(0, 10))

    match = set(range(10, 20))

    # load up the counter
    counter = CounterGather(query)
    with pytest.raises(ValueError):
        counter.add(match, "no name")

    assert counter.peek(query) == []


def test_counter_gather_empty_counter():
    # empty counter!
    query = set()
    counter = CounterGather(query)

    assert counter.peek(query) == []


def test_counter_gather_3_test_consume():
    # open-box testing of consume(...)
    query = set(range(0, 20))

    match1 = set(range(0, 10))
    match2 = set(range(7, 15))
    match3 = set(range(13, 17))


    # load up the counter
    counter = CounterGather(query)
    counter.add(match1, 'loc a')
    counter.add(match2, 'loc b')
    counter.add(match3, 'loc c')

    ### ok, dig into actual counts...
    import pprint
    pprint.pprint(counter.counter)
    pprint.pprint(counter.setlist)
    pprint.pprint(counter.locations)

    assert counter.setlist == [ match1, match2, match3 ]
    assert counter.locations == ['loc a', 'loc b', 'loc c']
    assert list(counter.counter.items()) == [(0, 10), (1, 8), (2, 4)]

    ## round 1

    cur_query = set(query)
    (cont, match, ident, intersect_mh) = counter.peek(cur_query)
    assert match == match1
    assert len(intersect_mh) == 10
    assert cur_query == query

    counter.consume(intersect_mh)
    assert counter.setlist == [ match1, match2, match3 ]
    assert counter.locations == ['loc a', 'loc b', 'loc c']
    assert list(counter.counter.items()) == [(1, 5), (2, 4)]

    ### round 2

    cur_query -= intersect_mh
    (cont, match, ident, intersect_mh) = counter.peek(cur_query)
    assert match == match2
    assert len(intersect_mh) == 5
    assert cur_query != query

    counter.consume(intersect_mh)
    assert counter.setlist == [ match1, match2, match3 ]
    assert counter.locations == ['loc a', 'loc b', 'loc c']
    assert list(counter.counter.items()) == [(2, 2)]

    ## round 3

    cur_query -= intersect_mh
    (cont, match, ident, intersect_mh) = counter.peek(cur_query)
    assert match == match3
    assert len(intersect_mh) == 2
    assert cur_query != query

    counter.consume(intersect_mh)
    assert counter.setlist == [ match1, match2, match3 ]
    assert counter.locations == ['loc a', 'loc b', 'loc c']
    assert list(counter.counter.items()) == []

    ## round 4 - nothing left!

    cur_query -= intersect_mh
    results = counter.peek(cur_query)
    assert not results

    counter.consume(intersect_mh)
    assert counter.setlist == [ match1, match2, match3 ]
    assert counter.locations == ['loc a', 'loc b', 'loc c']
    assert list(counter.counter.items()) == []
