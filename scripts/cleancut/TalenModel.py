### Copyright Daniel R. Chee 2017
### Altius Institute of Biomedical Sciences

import argparse
import collections
import intervaltree
import itertools
import numpy
import json
import re
import time
import logging

import bowtie
import genome
import primer
import file_parsers

log = logging.getLogger(__name__)

def talen_tree(chrom,start,end,seq_fasta=None,assembly="hg38",strand="+",ARRAY_MIN=15,ARRAY_MAX=21):
	"""
	"""
	log.debug(f"Making talen tree for {chrom}:{start}-{end} ({assembly} for lengths {ARRAY_MIN}-{ARRAY_MAX}), strand {strand}")
	if seq_fasta:
		fasta = file_parsers.load_fasta(seq_fasta)
		sequence = fasta[chrom][start:end]
	else:
		#sequence = genome.get_sequence(chrom,start,end,strand=strand,assembly=assembly)
		sequence = genome.get_sequence(chrom,start,end,assembly=assembly)
	### Defining REGEXs that identify talens on the positive and negative 
	### strands respectively.

	if strand == "+":
		tal_re = re.compile(f"[Tt]")
	else:
		tal_re = re.compile(f"[Aa]")

	seq_dict = {}
	coors = set()
	#### Finds all occurences of a talens across the sequence
	for tal in tal_re.finditer(sequence):
		i = tal.start() if strand == "+" else tal.start()-(ARRAY_MIN+1)
		j = tal.start()+(ARRAY_MIN+1) if strand == "+" else tal.start()+1
		while (i>=0) and (j<=len(sequence)) and ((j-i) < ARRAY_MAX+1):
			#tree[i:j] = sequence[i:j]
			seq_dict[(i,j)] = sequence[i:j]	
			coors.add((i,j))
			if strand == "+":
				j+=1
			else:
				i-=1

	numhits0 = bowtie.get_number_of_hits(seq_dict, assembly=assembly)
	
	tree = intervaltree.IntervalTree()
	for m0,m1 in coors:
		if strand == "+":
			tree[m0:m1] = TalenMonomer(chrom,start+m0,start+m1,sequence[m0+1:m1],sequence[m0],strand,numhits0[f"({m0}, {m1})"])
		else:
			tree[m0:m1] = TalenMonomer(chrom,start+m0,start+m1,sequence[m0:m1-1],sequence[m1-1],strand,numhits0[f"({m0}, {m1})"])

	return tree

def talen_dimer_iter(chrom,start,end,seq_fasta=None,assembly="hg38",strand="+",ARRAY_MIN=15,ARRAY_MAX=21,SPACER_MIN=14,SPACER_MAX=16):
	if seq_fasta:
		fasta = file_parsers.load_fasta(seq_fasta)
		sequence = fasta[chrom][start:end]
	else:
		sequence = genome.get_sequence(chrom,start,end,strand=strand,assembly=assembly)
		
	### Enumberating all talens for both strands in the region.
	pos_tree = talen_tree(chrom,start,end,seq_fasta=seq_fasta,assembly=assembly,strand="+",ARRAY_MIN=ARRAY_MIN,ARRAY_MAX=ARRAY_MAX)
	neg_tree = talen_tree(chrom,start,end,seq_fasta=seq_fasta,assembly=assembly,strand="-",ARRAY_MIN=ARRAY_MIN,ARRAY_MAX=ARRAY_MAX)
	### Building dimers from the extracted talens
	for l0,l1,left in sorted(pos_tree):
		for r0,r1,right in sorted(neg_tree[(l1+SPACER_MIN-1):(l1+SPACER_MAX+1)]):
			if (r0-l1) >= SPACER_MIN and (r0-l1) <= SPACER_MAX:
				# Correct the 
				yield TalenDimer(left, right, sequence[l1:r0], assembly)
	
def gc_content(sequence):
	"""
	"""
	score = sum([1 if base=="G" or base=="C" else 0 for base in sequence.upper()])/len(sequence)
	return score

def start_strength(sequence):
	"""
		Scores the start strength of the talen
		* i.e. the GC content of the first 5 bases
	"""
	score = sum([1 if base=="G" or base=="C" else 0 for base in sequence[:5].upper()])/len(sequence[:5])
	return score

def balance(sequence):
	"""
	"""
	x = numpy.array([n for (n,base) in enumerate(sequence) if (base=="G") or (base=="C")])
	if len(x) > 1:
		score = numpy.mean((x[1:]-x[:-1])-1)/len(sequence)
	else:
		score = 1.0
	return 1.0 - score

class TalenMonomer:
	
	def __init__(self, chromosome, start, end, sequence, flank, strand, numhits0, name="NA", rvds=None, RVD_map={"C":"HD","G":"NH","A":"NI","T":"NG"}, 
		backbone="pVAX",primer5=None,primer=None,score_GC=None,score_SS=None,score_balance=None,score_combined=None,
		numhits1=None,numhits2=None):

		self.name = name 
		self.chrom = chromosome
		self.i0 = start
		self.i1 = end
		self.sequence = sequence
		self.flank = flank 
		self.rvds = rvds if rvds else ",".join([RVD_map[base.upper()] for base in self.sequence])
		self.strand = strand # This could be either strand technically
		self.backbone = backbone
		self.primer = primer
		self.primer5 = primer5
		self.score_GC = score_GC if score_GC else gc_content(self.sequence)
		self.score_SS = score_SS if score_SS else start_strength(self.sequence)
		self.score_balance = score_balance if score_balance else balance(self.sequence)
		self.score_combined = score_combined if score_combined else numpy.mean([self.score_GC, self.score_SS, self.score_balance])
		self.numhits0 = numhits0
		self.numhits1 = None
		self.numhits2 = None
		return

	def json(self):
		return json.dumps(self.dict())

	def dict(self):
		return {
				"chromosome": self.chrom,
				"start": self.i0,
				"end": self.i1,
               	"sequence": self.sequence,
				"strand": self.strand,
				"name": self.name,
				"flank": self.flank,
				"primer": self.primer,
				"primer5": self.primer5,
				"rvds": self.rvds,
				"backbone": self.backbone,
				"score_GC": self.score_GC,
				"score_SS": self.score_SS,
				"score_balance": self.score_balance,
				"score_combined": self.score_combined,
				"numhits0": self.numhits0,
				"numhits1": self.numhits1,
				"numhits2": self.numhits2,
			}

class TalenDimer:
	
	def __init__(self,left,right,spacer,assembly,name="NA",RVDs={"C":"HD","G":"NH","A":"NI","T":"NG"},l_backbone="pVAX",r_backbone="pVAX"):
		# Talen monomers
		self.left = left
		self.right = right
		### Defining the information for the talen dimer
		self.dimer_name = name
		self.chrom = self.left.chrom
		self.spacer_0 = self.left.i1
		self.spacer_1 = self.right.i0
		self.score = (self.left.score_combined + self.right.score_combined)/2
		self.unique = int(f"{int(self.left.numhits0>1)}{int(self.right.numhits0>1)}", 2)
		self.assembly = assembly
		### TODO add assembly and strand
		self.relative_start = "NA" # we need to support the prefix ala +100, so a string for a number it is
		self.spacer = spacer.lower()
		self.cutsite = self.spacer_0 + ((self.spacer_1-self.spacer_0)//2)
		### Defining a unique value that can be used for hashing and comparing
		self.ident = (self.chrom, self.cutsite, f"T {self.left} {self.spacer} {self.right} A")

	def amplicon_sequence(self):
		if self.left.primer5 and self.right.primer5:
			return genome.get_sequence(self.chrom, self.left.primer5, self.right.primer5, \
				assembly=self.assembly, strand=self.left.strand)
		else:
			return genome.get_sequence(self.chrom, self.cutsite-100, self.cutsite+100, \
				assembly=self.assembly, strand=self.left.strand)
	
	def design_primers(self):
		if not(self.left.primer) or not(self.left.primer5):
			self.left.primer,self.left.primer5 = primer.left_closest(self.chrom, self.left.i0, assembly=self.assembly)
		if not(self.right.primer) or not(self.right.primer5):
			self.right.primer,self.right.primer5 = primer.right_farthest(self.chrom, self.cutsite+1, assembly=self.assembly)
		return

	def check_primers(self, num_cycles=150):
		if self.left.primer and self.left.primer5 and self.right.primer and self.right.primer5:
			return ((self.left.primer5 + len(self.left.primer) + num_cycles) - (self.right.primer5 + len(self.right.primer) - num_cycles)) > 10
		else:
			return False

	def unique_enough(self, min_hits=1, max_hits=10):
		return ((self.left.numhits0 <= min_hits) and (self.right.numhits0 <= max_hits)) or \
			((self.right.numhits0 <= min_hits) and (self.left.numhits0 <= max_hits))

	def set_relative_start(self, start):
		relative_start = self.cutsite-int(start)
		self.relative_start = f"+{relative_start}" if self.cutsite > int(start) else f"{relative_start}"

	def __eq__(self, other):
		return self.ident == other.ident
	
	def __hash__(self):
		return hash(self.ident)

	def __lt__(self,other):
		return self.ident < other.ident

	def simple_bed(self, strand="+"):
		"""
			Outputs the talen in 'simple' bed format with the genomic coordinates of the 
			cutsite and the sequence, either on the plus strand or minus strand, along with
			the respective number of hits for each talen.
		"""
		if strand == "+":
			return "\t".join([f"{self.chrom}", f"{self.cutsite}", f"{self.cutsite+1}", f"{self.dimer_name}", \
				f"T {self.left.sequence.upper()} {self.spacer.lower()} {self.right.sequence.upper()} A", \
				f"{self.left.numhits0}", f"{self.right.numhits0}"])
		else:
			return "\t".join([f"{self.chrom}", f"{self.cutsite}", f"{self.cutsite+1}", f"{self.dimer_name}", \
				(f"A {genome.complement(self.left.sequence.upper())} {genome.complement(self.spacer.lower())} " + \
				f"{genome.complement(self.right.sequence.upper())} T")[::-1], f"{self.left.numhits0}", f"{self.right.numhits0}"])

	def json(self):
		return json.dumps(self.dict())

	def dict(self):
		return {
			"dimer_name": self.dimer_name,
			"chromosome": self.chrom,
			"cutsite": self.cutsite,
			"relative_start": self.relative_start,
			"assembly": self.assembly,
			"spacer_sequence": self.spacer,
			"spacer_start": self.spacer_0,
			"spacer_end": self.spacer_1,
			"score": self.score,
			"left": self.left.dict(),
			"right": self.right.dict(),
		}

	def __str__(self):
		return "\t".join([f"{self.dimer_name}",f"{self.chrom}",f"{self.cutsite}",f"{self.cutsite+1}", \
			f"{self.assembly}",f"{self.spacer}",f"{self.score}",f"{self.unique}",f"{self.left.name}", \
			f"{self.chrom}",f"{self.left.i0}",f"{self.left.i1}",f"{self.left.strand}",f"{self.left.sequence}", \
			f"{self.left.flank}",f"{self.left.primer}",f"{self.left.primer5}",f"{self.left.rvds}", \
			f"{self.left.backbone}",f"{self.assembly}",f"{self.left.score_GC}",f"{self.left.score_SS}",\
			f"{self.left.score_balance}",f"{self.left.score_combined}",f"{self.left.numhits0}",\
			f"{self.left.numhits1}",f"{self.left.numhits2}",f"{self.right.name}", \
			f"{self.chrom}",f"{self.right.i0}",f"{self.right.i1}",f"{self.right.strand}",f"{self.right.sequence}", \
			f"{self.right.flank}",f"{self.right.primer}",f"{self.right.primer5}",f"{self.right.rvds}", \
			f"{self.right.backbone}",f"{self.assembly}",f"{self.right.score_GC}",f"{self.right.score_SS}",\
			f"{self.right.score_balance}",f"{self.right.score_combined}",f"{self.right.numhits0}",\
			f"{self.right.numhits1}",f"{self.right.numhits2}"])

def design_iter(dimer_path):
	"""
		Iterates over a tsv file of dimer designs that were output by in verbose form.
	"""
	for record in file_parsers.tsv_iter(dimer_path):
		### Creating a new talen dimer using info from the file
		left = TalenMonomer(record[1], int(record[10]), int(record[11]), record[13], record[14],record[12],int(record[24]))
		right = TalenMonomer(record[1],int(record[29]),int(record[30]),record[32],record[33],record[31],int(record[43]))
		dimer = TalenDimer(left,right,record[5],record[4])
		dimer.dimer_name = record[0]
		### Left enties that cannot be inferred when the object is created
		dimer.left.name = record[8]
		dimer.left.backbone = record[18]
		dimer.left.rvds = record[17]
		dimer.left.strand = record[12]
		dimer.left.primer = record[15] if record[15] != "None" else None 
		dimer.left.primer5 = int(record[16]) if record[16] != "None" else dimer.cutsite-100
		### Right enties that cannot be inferred when the object is created
		dimer.right.name = record[27]
		dimer.right.backbone = record[37]
		dimer.right.rvds = record[36]
		dimer.right.strand = record[31]
		dimer.left.primer = record[34] if record[34] != "None" else None 
		dimer.right.primer5 = int(record[35]) if record[35] != "None" else dimer.cutsite+100
		### Passing on the whole object
		yield dimer

def load_collection(talen_path):
	return {design.dimer_name:design for design in design_iter(talen_path)}

def build_verbose(sequence_path, assembly="hg38"):
	"""
		Currently this assumes that the all of the default stuff is the case
	"""
	### Empty dictionaries that will be passed to the bowtie module when filled.
	hits_dict = {}
	sequence_dict = {}
	### Filling the sequence dictionaries
	for name,sequence in file_parsers.tsv_iter(sequence_path):
		T,left,spacer,right,A = sequence.split(" ")
		sequence_dict[name] = sequence.replace(" ", "").upper()
		hits_dict[f"{name}_L"] = f"{T}{left}"
		hits_dict[f"{name}_R"] = f"{right}{A}"
	### Using bowtie to calculate the number of hits and coordinates
	coors = bowtie.get_coordinates(sequence_dict, assembly=assembly)
	numhits = bowtie.get_number_of_hits(hits_dict, assembly=assembly)
	### Re-iterating over the file to build the dimer objects
	for name,sequence in file_parsers.tsv_iter(sequence_path):
		T,left,spacer,right,A = sequence.split(" ")
		chrom,start,end,strand = coors[name]
		### All Talen objects will be output on the plus strand.
		plus_seq = genome.get_sequence(chrom,start,end,assembly=assembly)
		if strand == "+":
			l0,l1 = start,start+len(left)+1
			r0,r1 = start+len(left)+1+len(spacer),end
			dimer = TalenDimer(chrom,l0,l1,plus_seq[l0-start:l1-start],numhits[f"{name}_L"],
				r0,r1,plus_seq[r0-start:r1-start],numhits[f"{name}_R"], assembly)
		else:
			l0,l1 = start,start+len(right)+1
			r0,r1 = start+len(right)+1+len(spacer),end
			dimer = TalenDimer(chrom,l0,l1,plus_seq[l0-start:l1-start],numhits[f"{name}_R"],
				r0,r1,plus_seq[r0-start:r1-start],numhits[f"{name}_L"], assembly)
		### Filling in the names and printing in verbose format.
		dimer.dimer_name = name
		dimer.left_name = f"{name}_LEFT"
		dimer.right_name = f"{name}_RIGHT"
		print(dimer)
	return

def shared_amplicon(talens):
	dimers = sorted(list(talens))
	amplicon_left =  dimers[0].left_primer5 if dimers[0].left_primer5 else dimers[0].cutsite-150
	amplicon_right =  dimers[-1].right_primer5 if dimers[-1].right_primer5 else dimers[-1].cutsite+150
	amplicon = genome.get_sequence(dimers[0].chrom, amplicon_left, amplicon_right, \
		assembly=dimers[0].assembly, strand=dimers[0].left_strand)
	return dimers[0].chrom, amplicon_left, amplicon_right, amplicon

def setup_parser():

	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers()
	
	simple_output_parser = subparsers.add_parser("simple-output")
	simple_output_parser.add_argument("design_path")
	simple_output_parser.add_argument("--strand", type=str, default="+")
	simple_output_parser.set_defaults(which="simple-output")

	build_verbose_parser = subparsers.add_parser("build-verbose")
	build_verbose_parser.add_argument("sequence_path")
	build_verbose_parser.add_argument("--assembly", type=str, default="hg38")
	build_verbose_parser.set_defaults(which="build-verbose")

	check_primers_parser = subparsers.add_parser("check-primers")
	check_primers_parser.add_argument("design_path")
	check_primers_parser.add_argument("--strand", type=str, default="+")
	check_primers_parser.set_defaults(which="check-primers")

	output_primers_parser = subparsers.add_parser("output-primers")
	output_primers_parser.add_argument("design_path")
	output_primers_parser.set_defaults(which="output-primers")

	test_parser = subparsers.add_parser("test")
	test_parser.set_defaults(which="test")

	args = parser.parse_args()

	return args

def main(args):

	if args.which == "simple-output":
		for design in sorted(design_iter(args.design_path)):
			print(design.simple_bed(strand=args.strand))
	if args.which == "check-primers":
		for design in sorted(design_iter(args.design_path)):
			print(f"{design.dimer_name}: primers correct={design.check_primers()}")
	if args.which == "build-verbose":
		build_verbose(args.sequence_path, assembly=args.assembly)
	if args.which == "output-primers":
		illumina_left = "ACACGACGCTCTTCCGATCT"
		illumina_right = "GACGTGTGCTCTTCCGATCT"
		for design in sorted(design_iter(args.design_path)):
			print(f"{design.dimer_name}\t{illumina_left}NNNN{design.left_primer}\t{illumina_right}NNNN{design.right_primer}")

	return

if __name__ == "__main__":
	args = setup_parser()
	main(args)

