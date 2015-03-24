# This is a quick script to split up FASTQ files by barcode given
# Used to rescue tags from undeterminde state

import sys, os, gzip, re, operator
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import logging

import argparse


def parseArgs():
    parser = argparse.ArgumentParser(description='Split up fastq files by barcode')
    parser.add_argument('--mismatches', type=int, default=0, help='number of mismatches')
    parser.add_argument('--processing', dest='processing_file', action='store', required=True,
            help='processing.json to use (mandatory)')
    parser.add_argument('--suffix', dest='suffix', default='', help='suffix to add to sample names')
    parser.add_argument('--lane', dest='lane', type=int, default=1, help='Lane to process (default 1)')
    parser.add_argument('--autosuffix', action="store_true", default=False, help='Automatically guess a suffix name')
    parser.add_argument('infile', nargs='+')

    args = parser.parse_args()
    return args



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
    _ L \d{3}       # lane number again
    (                 # Start suffix
      _ R \d          # Read number
      _   \d{3}       # Count
    )                 # End suffix
    .fastq.gz
    $
    """, re.X)

    match = regex.search(filename)
    if match:
        return match.group(1)

    return ""

import json
def parse_processing_file(file, mismatches, suffix, lane):
    barcodes = dict()
    labels = dict()
    with open(file) as data_file:
        data = json.load(data_file)

    lane_libraries = [ l for l in data['libraries'] if l['lane'] == lane ]

    for library in lane_libraries:
        label = library['alignments'][0]['sample_name']
        barcode_indices = library['barcode_index'].split("-")
        barcode1 = barcode_indices[0]
        barcode2 = barcode_indices[1] if len(barcode_indices) > 1 else "" 

        lengths.add(( len(barcode1), len(barcode2) ))

        for b1 in mismatch(barcode1, mismatches):
            for b2 in mismatch(barcode2, mismatches):
                barcode = (b1, b2)
                # TODO: This can be smarter
                if barcode in barcodes:
                    print "Error! Barcode %s already taken, lower mismatches! (from %s+%s)" % (barcode, barcode1, barcode2)
                    sys.exit(1)
                barcodes[barcode] = label

        labels[label] = { "filtered": 0, "unfiltered": 0, "total": 0 }
        labels[label]["outfile"] = gzip.open("%s%s.fastq.gz" % (label, suffix), 'wb')

    return barcodes, labels

def split_file(filename, barcodes, labels):

    barcode_re = re.compile(r"""
            [012]:                  #
            ([YN]):                 # Fail/pass chastity filtering
            [01]:                   #
            ( [AGCTN] {6,8} )       # First barcode
            \+?                     # Optional separator (+)
            ( [AGCTN] {6,8} )?      # Optionally, second barcode
            $
            """, re.X)

    tally = 0
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
            print "Record %d in %s" % (tally, filename)
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

            # Replace recorded barcode
            sepIndex = record.rfind(':')
            record = record[:sepIndex + 1] + barcode1 + "+" + barcode2
            #write to FASTQ
            labels[label]['outfile'].write('@%s\n%s\n+\n%s\n' % (record, seq, qual))

            if filter == "Y":
                labels[label]["filtered"] += 1
            else:
                labels[label]['unfiltered'] += 1

    parsein.close()


def main(argv):
    args = parseArgs()

    global lengths
    lengths = set([])

    if args.autosuffix:
        args.suffix = guess_suffix(args.infile[0])
        print("Setting suffix to %s", args.suffix)

    barcodes, labels = parse_processing_file(args.processing_file, args.mismatches, args.suffix, args.lane)

    for filename in args.infile:
        split_file(filename, barcodes, labels) 

    print "BARCODE MATCHING"

    for label, info in labels.iteritems():
        print "%s\t%s" % (label, str(info))

    for label, info in labels.iteritems():
        info['outfile'].close()

if __name__ == "__main__":
    main(sys.argv)
