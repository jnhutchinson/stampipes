#!/usr/bin/env python

# Usage: $0 barcode1 barcode2 [...]
# Primarily intended for use in other scripts, primarly flowcells/setup.sh
# Feel free to make it read from a pipe as well.
# Algorithm copy&pasted from the $STAMPY version.

import os, sys, logging, re, math
from optparse import OptionParser
#from StamPy.samplesheet import SampleSheet
import itertools
import json

def hamming(str1, str2):
    """Returns the distance between two strings"""
    return sum(itertools.imap(unicode.__ne__, str1, str2))

def min_distance(barcodes):
    """Returns the min distance between a set of barcodes"""
    min_distance1 = 100
    min_distance2 = 100

    for index1 in range(0, len(barcodes)):
        for index2 in range(index1+1, len(barcodes)):
            distance1 = hamming(barcodes[index1][0], barcodes[index2][0])
            distance2 = hamming(barcodes[index1][1], barcodes[index2][1])
            #logging.debug("%s (%d) vs %s (%d): %d" % (barcodes[index1], index1, barcodes[index2], index2, distance))
            if distance1 < min_distance1:
                min_distance1 = distance1
            if distance2 < min_distance2:
                min_distance2 = distance2

    #logging.debug("Max distance: %d" % min_distance)
    print  (min_distance1, min_distance2)
    return (min_distance1, min_distance2)

def max_mismatch(barcodes):
    """Returns the max mismatch between a set of barcodes.  It is the max distance divided by two, then minus one."""

    distances = min_distance(barcodes)
    max_mismatch_distance = int(math.ceil(distances[0] / 2.0)) - 1
    logging.debug("Max mismatch: %d" % max_mismatch_distance)
    return max_mismatch_distance

def gen_snps(word, mismatches):
    for d in range(mismatches+1):
        for locs in itertools.combinations(range(len(word)), d):
            thisWord = [[char] for char in word]
            for loc in locs:
                origChar = word[loc]
                thisWord[loc] = [l for l in "ACGTN" if l != origChar]
            for poss in itertools.product(*thisWord):
                yield "".join(poss)

def generate_barcodes(barcode_tuple, mismatch_tuple):
    return set(itertools.product( gen_snps(barcode_tuple[0], mismatch_tuple[0]),
                              gen_snps(barcode_tuple[1], mismatch_tuple[1])))

def is_mismatch_level_okay(barcodes, mismatch_tuple):
    barcode_collection = set()
    for barcode in barcodes:
        new_barcodes = generate_barcodes( barcode, mismatch_tuple )
        #print barcode
        #print new_barcodes
        if barcode_collection.isdisjoint( new_barcodes ):
            barcode_collection.update(new_barcodes)
        else:
            return False
    return True

def get_max_mismatch_level(barcodes):
    mismatch_choices = [(1, 1), (1, 0), (0, 1), (0, 0)]
    for choice in mismatch_choices:
        no_collisions = is_mismatch_level_okay(barcodes, choice)
        if no_collisions:
            return choice
    return None

# Takes in string of mask & barcode, returns list of trimmed barcodes
# NB: This assumes that all index reads will start with an i, and be followed by one or more digits
#     e.g: i6n will work, but iiiiiin and i2n2i2 will not.
def apply_mask(mask, barcode_string):
    orig_barcodes = barcode_string.split('-')
    print(orig_barcodes)
    index_lengths = map(int, re.findall( r"""(?: i ( \d* ) ) """, mask, re.X | re.I))

    # Truncate barcodes to mask length
    #barcodes = [ b[:index_lengths[i]] for (i, b) in enumerate(orig_barcodes) ]
    barcodes = [ orig_barcodes[i][:l] for (i, l) in enumerate(index_lengths) ]
    return barcodes


def main(args = sys.argv):
    process_json = open('processing.json')
    data = json.load(process_json)
    process_json.close()

    mask = data['alignment_group']['bases_mask']

    lanes = {}
    final_mm = (100,100)
    for library in data['libraries']:
        lane = library['lane']
        barcodes = apply_mask(mask, library['barcode_index'])
        if len(barcodes) < 2:
            barcodes.append(u'')
        if lane not in lanes:
            lanes[lane] = set()
        lanes[lane].add(tuple(barcodes))
    for lane in lanes.keys():
        mm = get_max_mismatch_level(lanes[lane])
        print(lanes[lane])
        print "L%03d mm" % lane, mm
        final_mm = (min(mm[0], final_mm[0]), min(mm[1], final_mm[1]))

    print(",".join(map(str, final_mm)))

if __name__ == "__main__":
    main()
