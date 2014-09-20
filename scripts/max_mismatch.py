#!/usr/bin/env python

# Usage: $0 barcode1 barcode2 [...]
# Primarily intended for use in other scripts, primarly flowcells/setup.sh
# Feel free to make it read from a pipe as well.
# Algorithm copy&pasted from the $STAMPY version.

import os, sys, logging, re, math
from optparse import OptionParser
from StamPy.samplesheet import SampleSheet
import itertools

def hamming1(str1, str2):
    """Returns the distance between two strings"""
    return sum(itertools.imap(str.__ne__, str1, str2))

def max_distance(barcodes):
    """Returns the max distance between a set of barcodes"""
    max_distance = len(barcodes[0])

    for index1 in range(0, len(barcodes)):
        for index2 in range(index1+1, len(barcodes)):
            distance = hamming1(barcodes[index1], barcodes[index2])
            logging.debug("%s (%d) vs %s (%d): %d" % (barcodes[index1], index1, barcodes[index2], index2, distance))
            if distance < max_distance:
                max_distance = distance

    logging.debug("Max distance: %d" % max_distance)
    return max_distance

def max_mismatch(barcodes):
    """Returns the max mismatch between a set of barcodes.  It is the max distance divided by two, then minus one."""

    max_mismatch_distance = int(math.ceil(max_distance(barcodes) / 2.0)) - 1
    logging.debug("Max mismatch: %d" % max_mismatch_distance)
    return max_mismatch_distance

def main(args = sys.argv):
    print(max_mismatch(args))

if __name__ == "__main__":
    main()
