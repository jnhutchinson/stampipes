import sys
import argparse
import json

def parseArgs():
    parser = argparse.ArgumentParser(description='Split up fastq files by barcode')
    parser.add_argument('--processing', dest='processing_file', action='store', required=True,
            help='processing.json to use (mandatory)')
    parser.add_argument('--barcodes', dest='barcodes_file', action='store',required=True,
            help='barcode output to compare')
    args = parser.parse_args()
    return args

def main(argv):
    args = parseArgs()
    barcodes = json.load(open(args.barcodes_file))
    process = json.load(open(args.processing_file))

    success = "TRUE"

    for lib in process['libraries']:
        bc1 = (lib['barcode1']['reverse_sequence'])
        bc2 = (lib['barcode2']['sequence'])
        bcs = bc1 + bc2
        lane = lib['lane']

        for l in barcodes['Lanes']:
            if lane == l['LaneIndex']:
                if bcs in l['Counts']:
                    next
                else:
                    success = "FALSE"

    print(success)

if __name__ == "__main__":
    main(sys.argv)
