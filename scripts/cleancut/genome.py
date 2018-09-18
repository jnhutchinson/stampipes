### Copyright Daniel R. Chee 2017
### Altius Institute of Biomedical Sciences

import pyfaidx
import argparse
import os
import json
import logging

import file_parsers

log = logging.getLogger(__name__)

config_file = os.environ.get("TALEN_CONFIG", None)

def load_genomes(config_file):
	"""Load the genomes from the given config file, as {assembly: {"fasta": pyfaidx.Fasta(fasta_file), "faidx": pyfaidx.Faidx(fasta_file)}}"""
	with open(config_file) as raw_config:
		config = json.loads(raw_config.read())
		genome_fastqs = config.get("genome_fastq", {})
		return dict([(genome, {"fasta": pyfaidx.Fasta(fasta_file), "faidx": pyfaidx.Faidx(fasta_file)})
			for genome, fasta_file in genome_fastqs.items()])

	raise Exception ("Cannot load config file {}".format(config_file))

genome = load_genomes(config_file)

def complement(sequence):
	comp = {"A":"T", "T":"A", "G":"C", "C":"G","N":"N",
		"a":"t", "t":"a", "g":"c", "c":"g","n":"n"}
	return "".join(comp[base] for base in sequence)

def reverse_complement(sequence):
	comp = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N",
		"a":"t", "t":"a", "g":"c", "c":"g", "n":"n"}
	return "".join(comp[base] for base in sequence[::-1])

def get_sequence(chrom,start,end, assembly="hg38", strand="+"):
	"""From a given assembly (default hg38), retrieve the region from the genome of chrm:start-end, reverse complemented
	if strand is "-"."""

	if not assembly in genome:
		raise ValueError(f"Assembly {assembly} is not configured!")

	fasta = genome[assembly]["fasta"]
	faidx = genome[assembly]["faidx"]

	if not chrom in faidx.index:
		raise ValueError(f"Chromosome {chrom} is not in assembly {assembly}")

	if start < 0:
		log.warning(f"Talen tree received start under 0 for {chrom} ({assembly}): {start}, correcting to 0")
		start = 0

	chrom_length = faidx.index[chrom].rlen

	if end > chrom_length:
		log.warning(f"get_sequence received end over {chrom}'s length of {chrom_length}: {end}, setting to chrom length")
		end = chrom_length

	if strand == "+":
		return str(fasta[chrom][int(start):int(end)]).upper()
	else:
		return reverse_complement(str(fasta[chrom][int(start):int(end)])).upper()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers()

	get_sequence_parser = subparsers.add_parser("get-sequence")
	get_sequence_parser.add_argument("chrom")
	get_sequence_parser.add_argument("start")
	get_sequence_parser.add_argument("end")
	get_sequence_parser.add_argument("--strand", type=str, default="+")
	get_sequence_parser.add_argument("--assembly", type=str, default="hg38")
	get_sequence_parser.set_defaults(which="get-sequence")

	bed_to_fasta_parser = subparsers.add_parser("bed-to-fasta")
	bed_to_fasta_parser.add_argument("regions")
	bed_to_fasta_parser.set_defaults(which="bed-to-fasta")

	args = parser.parse_args()

	if args.which == "get-sequence":
		seq = get_sequence(args.chrom,args.start,args.end,assembly=args.assembly,
			strand=args.strand)
		print(f">{args.chrom}:{args.start}-{args.end}")
		print(seq)
	elif args.which == "bed-to-fasta":
		#for chrom,start,end,*misc in file_parser.bedN_iter()
		pass

