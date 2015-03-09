# This is a quick script to split up FASTQ files by barcode given
# Used to rescue tags from undeterminde state

# usage: splitbarcodes.py FILE1 FILE2 FILE3 FILE4 ...

import sys, os, gzip, re, operator
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import logging

import argparse
parser = argparse.ArgumentParser(description='Split up fastq files by barcode')
parser.add_argument('--mismatches', type=int, default=0, help='number of mismatches')
parser.add_argument('--barcodes', dest='barcodelist_file', action='store', required=True,
        help='barcode file (mandatory)')
parser.add_argument('--suffix', dest='suffix', default='', help='suffix to add to sample names')
parser.add_argument('--autosuffix', action="store_true", default=False, help='Automatically guess a suffix name')
parser.add_argument('infile', nargs='+')

args = parser.parse_args()

barcode_re = re.compile(r"""
        [012]:                  #
        ([YN]):                 # Fail/pass chastity filtering
        [01]:                   #
        ( [AGCTN] {6,8} )       # First barcode
        \+? ( [AGCTN] {6,8} )?  # Optionally, separator (+) and second barcode
        $
        """, re.X)

import itertools
# Generator for mismatches
# returns original string + whatever variations have at most [mismatches] hamming distance
def mismatch(word, mismatches):
    for d in range(mismatches+1):
        for locs in itertools.combinations(range(len(word)), d):
            thisWord = [[char] for char in word]
            for loc in locs:
                origChar = word[loc]
                thisWord[loc] = [l for l in "ACGTN" if l != origChar]
            for poss in itertools.product(*thisWord):
                yield "".join(poss)    

def guess_suffix(filename):
    regex = re.compile(r"""
    ^                 # Start
    (?: lane \d+ _ )? # optional lane number
    Undetermined      # TODO: Maybe change someday?
    (?: _S0 )?        # Optional index
    (                 # Start suffix
      _ L \d{3}       # lane number again
      _ R \d          # Read number
      _   \d{3}       # Count
    )                 # End suffix
    .fastq.gz
    $
    """, re.X)

    print filename
    match = regex.search(filename)
    if match:
        return match.group(1)

    return ""

global lengths
lengths = set([])

barcodes = dict()
labels = dict()

if args.autosuffix:
    args.suffix = guess_suffix(args.infile[0])
    print("Setting suffix to %s", args.suffix)

for line in open(args.barcodelist_file, 'r'):
    vals = line.strip().split("\t")
    label = vals[0]
    barcode1 = vals[1]
    if (len(vals) > 2):
        barcode2 = vals[2]
    else:
        if '-' in barcode1:
            b = barcode1.split('-')
            barcode1 = b[0]
            barcode2 = b[1]
        else:
            barcode2 = ""

    lengths.add(( len(barcode1), len(barcode2) ))

    # TODO: This doesn't work right yet! Only works for mm=0
    for b1 in mismatch(barcode1, args.mismatches):
        for b2 in mismatch(barcode2, args.mismatches):
            barcode = (b1, b2)
            # TODO: This can be smarter
            if barcode in barcodes:
                print "Error! Barcode %s already taken! (from %s+%s)" % (barcode, barcode1, barcode2)
                sys.exit(1)
            barcodes[barcode] = label

    labels[label] = { "filtered": 0, "unfiltered": 0, "total": 0 }
    labels[label]["outfile"] = gzip.open("%s%s.fastq.gz" % (label, args.suffix), 'wb')

tally = 0

def split_file(filename):
    global tally, labels, barcodes, barcode_re

    print "Splitting up file: %s" % filename

    if filename.endswith('.gz'):
        parsein = gzip.open(filename)
    else:
        parsein = open(filename, 'rU')

    for record, seq, qual in FastqGeneralIterator(parsein):
        tally += 1
        match = barcode_re.search(record)

        if not match:
            print "Could not match %s" % record
            print str(seq)
            print "Record %d" % (tally)
            sys.exit(1)

        matches = match.groups()
        filter = matches[0]
        barcode1 = matches[1]
        if (len(matches) > 2):
            barcode2 = matches[2]
        if barcode2 == None:
            barcode2 = ""
        matched=False
        for format in lengths:
            barcode = ( barcode1[:format[0]], barcode2[:format[1]] )
            if barcode in barcodes.keys():
                label = barcodes[barcode]
                matched=True
                break

        if matched:
            labels[label]['total'] += 1

            #write to FASTQ
            labels[label]['outfile'].write('@%s\n%s\n+\n%s\n' % (record, seq, qual))

            if filter == "Y":
                labels[label]["filtered"] += 1
            else:
                labels[label]['unfiltered'] += 1

    parsein.close()


[split_file(filename) for filename in args.infile]

print "BARCODE MATCHING"

for label, info in labels.iteritems():
    print "%s\t%s" % (label, str(info))

for label, info in labels.iteritems():
    info['outfile'].close()
