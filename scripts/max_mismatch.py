#!/usr/bin/env python

# Usage: $0
# Primarily intended for use in other scripts, primarly flowcells/setup.sh
# Algorithm design:
#   Mask barcodes according to use_bases_mask
#   Generate all allowable mismatch settings
#   For each setting (most permissive to least), generate all possible barcodes for each lane
#     If there are any collisions, check the next tighter mismatch setting

import os, sys, logging, re, math
from optparse import OptionParser
import itertools
import json

MAX_MISMATCH_LEVEL = 1  # Nextseq can allow 2, Hiseq 2500 allows only 1
POSSIBLE_MISMATCH_LEVELS = range( MAX_MISMATCH_LEVEL, -1, -1 )

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
    return set(itertools.product(* [
        gen_snps( barcode_tuple[i], mismatch_tuple[i] )
        for i in range(len(barcode_tuple)) ]))

def is_mismatch_level_okay(barcodes, mismatch_tuple):
    barcode_collection = set()
    for barcode in barcodes:
        new_barcodes = generate_barcodes( barcode, mismatch_tuple )
        if barcode_collection.isdisjoint( new_barcodes ):
            barcode_collection.update(new_barcodes)
        else:
            return False
    return True

def get_max_mismatch_level(lane_set, index_count):
    mismatch_choices = list( itertools.product( POSSIBLE_MISMATCH_LEVELS, repeat=index_count))
    lanes = lane_set.itervalues() # We don't actually care about lane labels (the key)
    for choice in mismatch_choices:
        no_collisions = all( [
            is_mismatch_level_okay(barcodes, choice)
            for barcodes in lanes ] )
        if no_collisions:
            return choice
    return None

# Takes in mask & barcode, returns list of trimmed barcodes
def apply_mask(mask, barcode_string):
    orig_barcodes = barcode_string.split('-')
    while len(orig_barcodes) < len(mask):
        orig_barcodes.append(u'')
    barcodes = [ orig_barcodes[i][:l] for (i, l) in enumerate(mask) ]
    return barcodes

# NB: This assumes that all index reads will start with an i, and be followed by one or more digits
#     e.g: i6n will work, but iiiiiin and i2n2i2 will not.
def parse_bases_mask(mask_string):
    mask = map(int, re.findall( r"""(?: i ( \d* ) )""", mask_string, re.X | re.I))
    return mask

def main(args = sys.argv):
    process_json = open('processing.json')
    data = json.load(process_json)
    process_json.close()

    mask_string = data['alignment_group']['bases_mask']
    mask = parse_bases_mask(mask_string)

    lanes = {}
    for library in data['libraries']:
        lane = library['lane']
        barcodes = apply_mask(mask, library['barcode_index'])
        if lane not in lanes:
            lanes[lane] = set()
        lanes[lane].add(tuple(barcodes))

    mismatch_level = get_max_mismatch_level( lanes, len(mask) )
    if not mismatch_level:
        sys.stderr.write("No allowable mismatch levels found, barcode collision?\n")
        sys.exit(1)
    print ",".join(map(str, mismatch_level))

if __name__ == "__main__":
    main()
