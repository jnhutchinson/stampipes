#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import itertools
import logging
import string
import sys

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

UMI_LEN = 8       # Fixed length for takara RNA UMTs
TRIM_PADDING = 4  # Arbitrary constant for trimming R1
TRIM_LEN = UMI_LEN + TRIM_PADDING


def parse_args():
    """ Just parse the args """
    parser = argparse.ArgumentParser(
        description='Annotate read names with Takara v3 UMT')
    parser.add_argument('r1_fastq')
    parser.add_argument('r2_fastq')
    parser.add_argument('out_r1')
    parser.add_argument('out_r2')

    args = parser.parse_args()
    return args


def attach_umt(r1, r2):
    """ Attach UMT to r1 & r2 and remove from sequences """
    # Put UMT in names
    umt = r2.seq[:UMI_LEN]
    umt_add = "#%s" % (umt)
    r1.id += umt_add
    r1.name = ""
    r1.description = " ".join(r1.description.split()[1:])
    r2.id += umt_add
    r2.description = " ".join(r2.description.split()[1:])
    r2.name = ""

    # Check for presence of UMT in R1 - this indicates a short
    # fragment that needs to be trimmed

    # Slide the end of r2's revcom around until it aligns
    # e.g: i = 2:
    #                     |--needle---|
    # r1     AAAAAAAAAAAAA AAGG UMTUMT
    #           read      |pad-|r1-UMT|
    # r2_rev NNNAAAAAAAAAA AAGG UMTUMT UM
    #                     |--haystack-|
    #                          |-r2-UMT-|
    needle = str(r1.seq)[-TRIM_LEN:]
    haystack = str(r2.seq[:TRIM_LEN+UMI_LEN].reverse_complement())
    pos = haystack.find(needle)
    if pos > 0:
        r1 = r1[:-pos]

    # Remove from r2
    r2 = r2[UMI_LEN:]
    return (r1, r2)


def main(argv):
    args = parse_args()
    logging.basicConfig(level=logging.WARN, format=log_format)

    with open(args.r1_fastq) as r1_in, \
            open(args.r2_fastq) as r2_in, \
            open(args.out_r1, 'wt') as r1_out, \
            open(args.out_r2, 'wt') as r2_out:

        r1_seqIO = SeqIO.parse(r1_in, "fastq")
        r2_seqIO = SeqIO.parse(r2_in, "fastq")
        try:
            while True:
                (r1, r2) = attach_umt(next(r1_seqIO), next(r2_seqIO))
                r1_out.write(r1.format("fastq"))
                r2_out.write(r2.format("fastq"))
        except StopIteration:
            logging.debug("EOF reached")

if __name__ == "__main__":
    main(sys.argv)
