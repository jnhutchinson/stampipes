#!/usr/bin/env python

# TODO: Doesn't work quite yet!!!!!

import argparse
import sys
import pysam
from collections import defaultdict

# Some bad globals
read_histo = defaultdict(int)
pair_map = {}


def parser_setup():

    script_options = {
        "infile": "/dev/stdin",
        "outfile": "/dev/stdout",
        "histfile": None,
    }

    parser = argparse.ArgumentParser()

    parser.add_argument("--hist", dest="histfile",
                        help="Write histogram of duplicates to this file")

    parser.add_argument("-i", "--infile", dest="infile",
                        help="Read from this file")
    parser.add_argument("-o", "--outfile", dest="outfile",
                        help="Write to this file")

    parser.set_defaults(**script_options)
    parser.set_defaults(quiet=False, debug=False)

    return parser


# Takes a list of reads, returns one sorted such that the "best" read is at the
# top.  Currently, that means highest mapping quality, with ties broken by
# query_name
def read_sort_key(x):
    return (-1 * x.mapping_quality, x.query_name)


# Takes in a list of reads, returns a list of lists - one list for each unique
# read. These reads should have identical chromosome & reference-start fields.
# They are sorted so that the highest-mapping quality read is first in each
# sub-list.
def find_dups(reads):
    # Don't waste time if there's only one read at this position
    if len(reads) == 1:
        return [reads]

    lists = defaultdict(list)
    lists2 = defaultdict(list)
    for r in reads:
        try:
            key = (r.template_length, r.get_tag("XD:Z"))
        except KeyError:
            key = r.template_length
        lists[key].append(r)

    return [sorted(l, key=read_sort_key) for l in lists.values()]


# Sets a read's duplicate flag, returns it
def set_dup(read, is_dup):
    if is_dup:
        read.flag |= 1024
    else:
        read.flag &= ~1024
    return read


def process_chunk(chunk, output):
    new_reads = []
    already_seen_reads = []
    for r in chunk:
        if r.query_name in pair_map:
            already_seen_reads.append(r)
        else:
            pair_map[r.query_name] = None
            new_reads.append(r)

    read_sets = find_dups(new_reads)

    for l in read_sets:
        # Add to the histogram
        read_histo[len(l)] += 1
        # Mark all but the highest-quality read as duplicates
        l[0] = set_dup(l[0], False)
        pair_map[l[0].query_name] = False
        for r in l[1:]:
            r = set_dup(r, True)
            pair_map[r.query_name] = True

        # Print the read set
        for r in l:
            output.write(r)

    # Print out already seen reads
    for r in already_seen_reads:
        r = set_dup(r, pair_map[r.query_name])
        del pair_map[r.query_name]
        output.write(r)


def main(args=sys.argv):
    parser = parser_setup()
    poptions = parser.parse_args()

    input = pysam.AlignmentFile(poptions.infile, 'r')
    output = pysam.AlignmentFile(poptions.outfile, 'wb0', template=input)

    chr = None
    pos = None
    current_chunk = []
    for read in input:
        if read.reference_id != chr or read.reference_start != pos:
            process_chunk(current_chunk, output)
            current_chunk = []
            chr = read.reference_id
            pos = read.reference_start
        current_chunk.append(read)
    process_chunk(current_chunk, output)

    if poptions.histfile:
        with open(poptions.histfile, 'w') as h:
            for k in sorted(read_histo.keys()):
                h.write("%s\t%s\n" % (k, read_histo[k]))


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell
# after importing without automatically running it
if __name__ == "__main__":
    main()
