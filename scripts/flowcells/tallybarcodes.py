# usage: tallybarcodes.py FILE1 FILE2 FILE3 FILE4 ...

import sys, os, gzip, re, operator, argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def parser_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--passing-only", dest="passing_only", action="store_true",
        help="Only tally barcodes that passed chastity filtering.")
    parser.add_argument('infiles', metavar='N', nargs='+',
        help="Input fastq files to tally.")

    return parser

parser = parser_setup()
args = parser.parse_args()

tally = 0
filtered_tally = 0

barcodes = dict()
barcode_re = re.compile("[12]:([YN]):[01]:([AGCTN]+)$")

def tally_file(filename):
    global tally, filtered_tally, barcodes, barcode_re
    
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

        filter = match.groups()[0]
        barcode = match.groups()[1]

        if filter == "Y":
            filtered_tally += 1
            if args.passing_only:
                continue

        if not barcode in barcodes:
            barcodes[barcode] = 0
        barcodes[barcode] += 1

    parsein.close()

[tally_file(filename) for filename in args.infiles]

print "Total records: " + str(tally)
print "Failing filter: " + str(filtered_tally)

print "BARCODES"

for barcode, tally in sorted(barcodes.iteritems(), key=operator.itemgetter(1), reverse = True):
    print "%s\t%d" % (barcode, tally)
