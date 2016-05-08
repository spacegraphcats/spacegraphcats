#! /usr/bin/env python3
import screed
import argparse
import json

DEFAULT_KSIZE=32

def kmers(sequence, K):
    for i in range(len(sequence) - K + 1):
        yield sequence[i:i+K]

def linear_seg_yield(sequence, K, kmer_counts):
    kmer_list = list(kmers(sequence, K))
    ki = 0

    linear = kmer_list[0][:-1]
    while ki < len(kmer_list):
        kmer = kmer_list[ki]

        if kmer_counts[kmer] > 1:
            if len(linear) >= K:
                yield linear
            yield kmer
            linear = kmer[1:]
        else:
            linear += kmer[-1]

        ki += 1
    if linear:
        if len(linear) >= K:
            yield linear

def test_linear_seg_yield():
    testseq = "AAAAATTTTTCCCCCGGGGG"
    K=5

    kmer_counts = {}
    for kmer in kmers(testseq, K):
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1

    x = list(linear_seg_yield(testseq, 5, kmer_counts))
    assert x[0] == testseq

    kmer_counts["TTTTT"] = 2
    x = list(linear_seg_yield(testseq, 5, kmer_counts))
    print(x)
    assert x == ['AAAAATTTT', 'TTTTT', 'TTTTCCCCCGGGGG']

    kmer_counts["TTTTT"] = 1
    kmer_counts["AAAAA"] = 2
    x = list(linear_seg_yield(testseq, 5, kmer_counts))
    assert x == ['AAAAA', 'AAAATTTTTCCCCCGGGGG']
    print(x)

    kmer_counts["GGGGG"] = 2
    kmer_counts["AAAAA"] = 1
    x = list(linear_seg_yield(testseq, 5, kmer_counts))
    print(x)
    assert x == ['AAAAATTTTTCCCCCGGGG', 'GGGGG']


def find_subseq(subseq, longer):
    loc = -1
    posns = []
    while 1:
        loc = longer.find(subseq, loc + 1)
        if loc == -1:
            break
        else:
            posns.append(loc)
    return posns

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqfiles', nargs='+')
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('-k', '--ksize', default=DEFAULT_KSIZE, type=int)
    parser.add_argument('--repeat-table', default=None)
    args = parser.parse_args()

    kmer_counts = {}
    linear_segments = {}
    linear_segments_r = {}
    linear_segments_origins = {}
    linear_segment_count = 1
    prev_dict = {}
    post_dict = {}
    is_repeat_dict = {}
    segment_order = {}

    assert len(args.seqfiles) == 1
    num_sequences = 0

    # first pass: compute k-mer counts
    print('loading from files', args.seqfiles)
    for filename in args.seqfiles:
        print('loading', filename)
        for record in screed.open(filename):
            num_sequences += 1
            for kmer in kmers(record.sequence, args.ksize):
                kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1

    assert num_sequences == 1
    #return kmer_counts

    if args.repeat_table:
        fp = open(args.repeat_table)
        for line in fp:
            start, end, feature = line.strip().split()
            assert feature == 'repeat_region'
            start = int(start)
            end = int(end)
        for kmer in kmers(record.sequence, args.ksize):
            is_repeat_dict[kmer] = True

    # second pass: compute linear segments
    pos = 0
    for filename in args.seqfiles:
        print('loading x 2', filename)

        for record in screed.open(filename):
            for linear in linear_seg_yield(record.sequence, args.ksize,
                                           kmer_counts):
                pos += 1
                if linear in linear_segments_r:
                    pass
                else:
                    this_line = linear_segment_count
                    linear_segments[this_line] = linear
                    linear_segments_r[linear] = this_line
                    linear_segments_origins[this_line] = record.name

                    linear_segment_count += 1

                this_line = linear_segments_r[linear]
                x = segment_order.get(this_line, [])
                x.append(pos)
                segment_order[this_line] = x

    # third pass: compute graph from linear segments & k-mer connections
    for filename in args.seqfiles:
        print('loading x 3', filename)

        for record in screed.open(filename):
            segs = linear_seg_yield(record.sequence, args.ksize, kmer_counts)
            segs = list(segs)
            print(record.name, len(segs))

            for si in range(len(segs)):
                this_seg = segs[si]
                this_id = linear_segments_r[this_seg]
                if si > 0:
                    prev = segs[si - 1]
                    prev_id = linear_segments_r[prev]

                    x = prev_dict.get(this_id, [])
                    if prev_id not in x:
                        x.append(prev_id)
                        prev_dict[this_id] = x

                if si < len(segs) -1:
                    post = segs[si + 1]
                    post_id = linear_segments_r[post]

                    x = post_dict.get(this_id, [])
                    if post_id not in x:
                        x.append(post_id)
                        post_dict[this_id] = x

    # convert to JSON
    if args.output:
        fp = open(args.output, 'w')
        output = []
        for linear_id in linear_segments:
            d = {}
            d['segment_id'] = linear_id
            d['end1'] = prev_dict.get(linear_id, [])
            d['end2'] = post_dict.get(linear_id, [])
            d['segment_size'] = len(linear_segments[linear_id])
            d['segment_order'] = segment_order[linear_id]
            d['multiplicity'] = 1
            d['is_repeat'] = 0
            if len(linear_segments[linear_id]) == args.ksize:
                kmer = linear_segments[linear_id]
                d['multiplicity'] = kmer_counts[kmer]
                d['is_repeat'] = is_repeat_dict.get(kmer, 0)

            if 0 and linear_id == 878 or linear_id == 905:
                print('YO!!!')
                print(d)

            output.append(d)

        fp.write(json.dumps(output))
        fp.close()

    print('computed %d segments' % len(linear_segments))

    mult1_total = 0
    multp_total = 0
    for linear_id in linear_segments:
        assert len(linear_segments[linear_id]) >= args.ksize
        prev_l = prev_dict.get(linear_id, [])
        post_l = post_dict.get(linear_id, [])
        x = set(prev_l)
        x.update(post_l)

        if len(linear_segments[linear_id]) > args.ksize:
            mult1_total += len(linear_segments[linear_id])
        else:
            multp_total += 1

    print('multiplicity one chains:', mult1_total)
    print('high multiplicity k-mers:', multp_total)
    print('total sequence length:', len(record.sequence))


    if 0:
        seq = linear_segments[71]
        print('71 is:', seq)
        print(find_subseq(seq, record.sequence))

        seq = linear_segments[881]
        print('881 is:', seq)
        print(find_subseq(seq, record.sequence))

        seq = linear_segments[878]
        print('878 is:', seq)
        print(find_subseq(seq, record.sequence))

        seq = linear_segments[905]
        print('905 is:', seq)
        print(find_subseq(seq, record.sequence))

if __name__ == '__main__':
    kmer_counts =  main()
