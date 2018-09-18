### Copyright Daniel R. Chee 2018
### Altius Institute of Biomedical Sciences

import argparse
import collections
import gzip
import itertools
import os
import regex

import file_parsers
import genome
import illumina

from Bio import SeqIO

index_re = regex.compile(r"(?P<i7>[ATCG]+)\+(?P<i5>[ATCG]+)")
class FastqEntry:
	
	def __init__(self,header,sequence,desc,score):
		self.header = header
		self.sequence = sequence
		self.desc = desc
		self.score = score
		self.index1,self.index2 = None,None
		for m in index_re.finditer(self.header):
			self.index1 = m.group("i7")
			self.index2 = m.group("i5")
		return

	def __str__(self):
		return f"{self.header}\n{self.sequence}\n{self.desc}\n{self.score}"

def fq_iter(fq_path):
	if fq_path.endswith(".gz"):
		fq_file = gzip.open(fq_path, "rt")
	else:
		fq_file = open(fq_path, "r")
	for n,line in zip(itertools.cycle(range(4)), fq_file):
		if n == 0:
			header=line.strip()
		elif n == 1:
			sequence = line.strip()
		elif n == 2:
			desc = line.strip()
		elif n == 3:
			score = line.strip()
			yield FastqEntry(header,sequence,desc,score)
	fq_file.close()

def demultiplex(forward_fq, reverse_fq, spec_sheet, outdir):
	"""
		Currently only looks for exact matches to the i7 and i5 indeces, theoretically 
		this could be retooled to allow for one mismatch, but that shouldn't really matter all
		that much. 
	"""
	### Creating the outdir if it doesn't exist
	if not os.path.isdir(outdir):
		os.makedirs(outdir)
	### Creating a dictionary of files
	forward_files = {}
	reverse_files = {}
	read_count = {}
	sample = {}
	for name,i7,i5,*etc in file_parsers.tsv_iter(spec_sheet):
		index_combo = (illumina.i7[i7], illumina.i5[i5])
		forward_files[index_combo] = open(f"{outdir}/{name}_r1.fastq", "w")
		reverse_files[index_combo] = open(f"{outdir}/{name}_r2.fastq", "w")
		read_count[index_combo] = 0
		sample[index_combo] = name
	### If it's an index combo we care about, put it in the appropriate file
	for read1,read2 in zip(fq_iter(forward_fq),fq_iter(reverse_fq)):
		index_combo = (read1.index1,read1.index2)
		if index_combo in forward_files:
			print(read1, file=forward_files[index_combo])
			print(read2, file=reverse_files[index_combo])
			read_count[index_combo]+=1
	### Close all of the files and write summary file
	with open(f"{outdir}/summary.txt", "w") as summary_file:
		print(f"name\ti7\ti5\t#reads", file=summary_file)
		for index_combo in forward_files: 
			forward_files[index_combo].close()
			reverse_files[index_combo].close()
			print(sample[index_combo], index_combo[0], index_combo[1], f"{read_count[index_combo]:,}", sep="\t", file=summary_file)
	return

def merge_paired_ends(forward_fq, reverse_fq):
	
	return

def combine_reads(merged_fq, reads_fa, unmerged_forward_fq=None, unmerged_reverse_fq=None):
	with open(reads_fa, "w") as reads_file:
		### merged fastq
		found = set()
		sequence_counts = collections.Counter()
		for record in itertools.islice(SeqIO.parse(merged_fq, "fastq"), None):
			umi = (str(record.seq[:4]), str(record.seq[-4:]))
			#if umi not in found:
			sequence_counts[str(record.seq[4:-4])]+=1
			found.add(umi)
		
		if unmerged_forward_fq:
			for record in itertools.islice(SeqIO.parse(unmerged_forward_fq, "fastq"), None):
				umi = (str(record.seq[:4]), "NNNN")
				#if umi not in found:
				sequence_counts[str(record.seq[4:])]+=1
				found.add(umi)

		if unmerged_reverse_fq:
			for record in itertools.islice(SeqIO.parse(unmerged_reverse_fq, "fastq"), None):
				umi = ("NNNN", str(record.seq[-4:]))
				#if umi not in found:
				sequence_counts[str(record.seq[4:])]+=1
				found.add(umi)

		for n,(seq,count) in enumerate(sequence_counts.most_common()):
			print(f">read_{n+1}_{count}\n{seq}", file=reads_file)
	return

def merge_selex(forward_fq, reverse_fq, out_path):
	selex_re = regex.compile("([ATCG]{18,25})AGATCGGAAGAGC")
	sequences = collections.Counter()
	i,j = 0,0
	for read1,read2 in zip(fq_iter(forward_fq),fq_iter(reverse_fq)):
		s1,s2 = "",""
		for m in selex_re.findall(read1.sequence):
			s1 = m
			i+=1
		for m in selex_re.findall(read2.sequence):
			s2 = m
			j+=1
		if s1:
			sequences[s1]+=1
		elif s2:
			sequences[s2]+=1
	seq_length = collections.Counter()
	with open(out_path, "w") as out_fa:
		for n,(seq,count) in enumerate(sequences.most_common()):
			print(f">read{n+1}_{count}\n{seq}", file=out_fa)
			seq_length[len(seq)]+=1
	print(f"{i:,}\t{j:,}\t{seq_length[18]:,}\t{seq_length[21]:,}\t{seq_length[25]:,}\t{sum(seq_length.values()):,}")
	return

def setup_parser():
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers(dest="command")
	###
	demultiplex_parser = subparsers.add_parser("demultiplex")
	demultiplex_parser.add_argument("forward_fq")
	demultiplex_parser.add_argument("reverse_fq")
	demultiplex_parser.add_argument("spec_sheet")
	demultiplex_parser.add_argument("outdir")
	###
	merged_paired_ends_parser = subparsers.add_parser("merge-paired-ends")
	merged_paired_ends_parser.add_argument("forward_fq")
	merged_paired_ends_parser.add_argument("reverse_fq")
	###
	combine_reads_parser = subparsers.add_parser("combine-reads")
	combine_reads_parser.add_argument("merged_fq")
	combine_reads_parser.add_argument("reads_fa")
	combine_reads_parser.add_argument("--unmerged-forward-fq", type=str, default=None)
	combine_reads_parser.add_argument("--unmerged-reverse-fq", type=str, default=None)
	###
	merge_selex_parser = subparsers.add_parser("merge-selex")
	merge_selex_parser.add_argument("forward_fq")
	merge_selex_parser.add_argument("reverse_fq")
	merge_selex_parser.add_argument("out_fa")
	###
	args = parser.parse_args()
	if args.command == "demultiplex":
		demultiplex(args.forward_fq, args.reverse_fq, args.spec_sheet, args.outdir)
	if args.command == "combine-reads":
		combine_reads(args.merged_fq, args.reads_fa, args.unmerged_forward_fq, args.unmerged_reverse_fq)
	if args.command == "merge-selex":
		merge_selex(args.forward_fq, args.reverse_fq, args.out_fa)
	return

if __name__ == "__main__":
	setup_parser()

