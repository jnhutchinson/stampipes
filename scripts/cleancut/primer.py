### Copyright Daniel R. Chee 2017
### Altius Institute of Biomedical Sciences

import re
import heapq
import primer3
import argparse
import itertools

import file_parsers
import bowtie
import genome

def tm(seq):
	return primer3.calcTm(seq, tm_method='breslauer', salt_corrections_method='schildkraut')

def left_closest(chrom, coor, min_tm=57, max_tm=63, min_len=15, max_len=25, min_dist=30, max_dist=100, assembly="hg38"):
	seq = genome.get_sequence(chrom,coor-(max_dist+max_len),coor-min_dist,assembly=assembly)
	primer_dict = {}
	for i in range(len(seq)):
		for j in range(i+min_len,i+max_len+1):
			if (len(seq[i:j]) >= min_len) and (len(seq[i:j]) <= max_len):
				cur_tm = tm(seq[i:j])
				if (cur_tm >= min_tm) and (cur_tm <= max_tm):
					primer_dict[f">primer_{i}_{j-i}"] = seq[i:j]
	if len(primer_dict) > 0:
		numhits = bowtie.get_number_of_hits(primer_dict)
		[primer_id] = heapq.nsmallest(1,primer_dict,key=lambda x: (numhits[x], abs((int(x.split("_")[1])+len(primer_dict[x]))-coor)))
		return primer_dict[primer_id],coor-(max_dist+max_len)+int(primer_id.split("_")[1])
	else:
		return None,None

def right_farthest(chrom, coor, min_tm=57, max_tm=63, min_len=15, max_len=25, min_dist=30, max_dist=100, assembly="hg38"):
	seq = genome.get_sequence(chrom,coor+min_dist,coor+(max_dist+max_len),assembly=assembly,strand="-")
	primer_dict = {}
	for i in range(len(seq)):
		for j in range(i+min_len,i+max_len+1):
			if (len(seq[i:j]) >= min_len) and (len(seq[i:j]) <= max_len):
				cur_tm = tm(seq[i:j][::-1])
				if (cur_tm >= min_tm) and (cur_tm <= max_tm):
					primer_dict[f">primer_{i}_{j-i}"] = seq[i:j]
	if len(primer_dict) > 0:
		numhits = bowtie.get_number_of_hits(primer_dict)
		[primer_id] = heapq.nlargest(1,primer_dict,key=lambda x: (numhits[x], abs((int(x.split("_")[1])-len(primer_dict[x]))-coor)))
		return primer_dict[primer_id],coor+min_dist+int(primer_id.split("_")[1])
	else:
		return None,None

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers()
	
	left_closest_parser = subparsers.add_parser("left-closest")
	left_closest_parser.add_argument("chromosome")
	left_closest_parser.add_argument("coordinate", type=int)
	left_closest_parser.add_argument("--min-tm", type=int, default=57)
	left_closest_parser.add_argument("--max-tm", type=int, default=63)
	left_closest_parser.add_argument("--min-len", type=int, default=15)
	left_closest_parser.add_argument("--max-len", type=int, default=25)
	left_closest_parser.add_argument("--min-dist", type=int, default=35)
	left_closest_parser.add_argument("--max-dist", type=int, default=100)
	left_closest_parser.add_argument("--assembly", type=str, default="hg38")
	left_closest_parser.set_defaults(which="left-closest")
		
	right_closest_parser = subparsers.add_parser("right-closest")
	right_closest_parser.add_argument("chromosome")
	right_closest_parser.add_argument("coordinate", type=int)
	right_closest_parser.add_argument("--min-tm", type=int, default=57)
	right_closest_parser.add_argument("--max-tm", type=int, default=63)
	right_closest_parser.add_argument("--min-len", type=int, default=15)
	right_closest_parser.add_argument("--max-len", type=int, default=25)
	right_closest_parser.add_argument("--min-dist", type=int, default=35)
	right_closest_parser.add_argument("--max-dist", type=int, default=100)
	right_closest_parser.add_argument("--assembly", type=str, default="hg38")
	right_closest_parser.set_defaults(which="right-closest")
		
	args = parser.parse_args()

	if args.which == "left-closest":
		primer = left_closest(args.chromosome,args.coordinate, min_tm=args.min_tm, max_tm=args.max_tm, \
						min_len=args.min_len, max_len=args.max_len, assembly=args.assembly)
	elif args.which == "right-closest":
		primer = right_closest(args.chromosome,args.coordinate, min_tm=args.min_tm, max_tm=args.max_tm, \
						min_len=args.min_len, max_len=args.max_len, assembly=args.assembly)

