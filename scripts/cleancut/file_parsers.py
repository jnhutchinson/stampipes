### Copyright Daniel R. Chee 2017
### Altius Institute of Biomedical Sciences

import csv
import gzip
import heapq
import itertools
import regex
import yaml

def csv_iter(csv_path, l0=0, l1=None):
	"""
		Generator wrapper for the csv parser. Punts opening and closing the file
		to the context manager.
	"""
	with open(csv_path, "rU") as csv_file:
		lines = (line.strip().split(",") for line in csv_file if line)
		yield from itertools.islice(lines, l0, l1)

def tsv_iter(tsv_path, l0=0, l1=None):
	"""
		Generator wrapper for a tsv parser. Punts opening and closing the file
		to the context manager.
	"""	
	try:
		with open(tsv_path, "rU") as tsv_file:
			lines = (line for line in csv.reader(tsv_file,delimiter="\t") if line)
			yield from itertools.islice(lines, l0, l1)
	except TypeError:
		lines = (line for line in csv.reader(tsv_path,delimiter="\t") if line)
		yield from itertools.islice(lines, l0, l1)
		

def fastq_iter(fastq_path):
	"""
		Iterator that yields just the sequence from a fastq file record by record.
	"""
	with open(fastq_path, 'rU') as fastq_file:
		lines = (line.strip() for line in fastq_file)
		yield from (line for i,line in enumerate(lines) if i%4 == 1)

count_re = regex.compile(r"read_[0-9]*_([0-9]*)")

class FastaEntry:
	
	def __init__(self, _ID, sequence):
		self.ID = _ID
		self.sequence = sequence
		return
	
	def __str__(self):
		return f">{self.ID}\n{self.sequence}"

def fa_iter(fa_path):
	with open(fa_path, "r") as fa_file:
		sequence = ""
		for line in fa_file:
			line = line.strip()
			if line.startswith(">"):
				if sequence:
					yield FastaEntry(_ID, sequence)
					sequence = ""
				_ID = line.replace(">", "")
			else:
				sequence = f"{sequence}{line}"
		if sequence:
			yield FastaEntry(_ID, sequence)

class FastqEntry:
	
	def __init__(self, meta, sequence, strand, score):
		self.meta = meta
		self.sequence = sequence
		self.strand = strand
		self.score = score
		return
	
	def __str__(self):
		return "{0}\n{1}\n{2}\n{3}".format(self.meta,self.sequence,self.strand,self.score)

def fq_iter(fastq_path, compressed=True):
	if compressed:
		with gzip.open(fastq_path, "rt") as fastq_file:
			meta,sequence,strand,score = "","","",""
			for n,line in enumerate(fastq_file):
				if n%4 == 0:
					meta = line.strip()
				elif n%4 == 1:
					sequence =line.strip()
				elif n%4 == 2:
					strand = line.strip()
				elif n%4 == 3:
					score = line.strip()
					yield FastqEntry(meta,sequence,strand,score)
	else:
		with open(fastq_path, "r") as fastq_file:
			meta,sequence,strand,score = "","","",""
			for n,line in enumerate(fastq_file):
				if n%4 == 0:
					meta = line.strip()
				elif n%4 == 1:
					sequence =line.strip()
				elif n%4 == 2:
					strand = line.strip()
				elif n%4 == 3:
					score = line.strip()
					yield FastqEntry(meta,sequence,strand,score)

class PslEntry:

	#def __init__(self, blat_record, ref_seq, query_seq):
	def __init__(self, blat_record):
		"""
			Python object representation of a BLAT alignment.
			Mostly to make it easier on the eyes to parse alignment.
		"""
		### Pulling all of the information from the
		self.matches = int(blat_record[0])
		self.mismatches = int(blat_record[1])
		self.num_inserts = int(blat_record[6])
		self.num_bases_inserted = int(blat_record[7])
		self.strand = blat_record[8]
		self.q_name = blat_record[9]
		self.r_name = blat_record[13]
		self.block_sizes = blat_record[18].split(',')[:-1]
		self.q_starts = blat_record[19].split(',')[:-1]
		self.r_starts = blat_record[20].split(',')[:-1]
		self.blat_record = blat_record
		self.q0 = int(self.q_starts[0])
		self.q1 = int(self.q_starts[-1]) + int(self.block_sizes[-1])
		self.r0 = int(self.r_starts[0])
		self.r1 = int(self.r_starts[-1]) + int(self.block_sizes[-1])
		self.count = int(count_re.search(blat_record[9]).group(1)) if count_re.match(blat_record[9]) else 1
		#self.score = (self.matches-self.mismatches)/(self.q1-self.q0)
		#self.score = (self.matches)/(self.q1-self.q0)
		### Setting the string equal to None
		self.string = "\t".join(blat_record)

	def get_alignment_groups(self):
		"""
			Generator function that iterates sequentially (from 'left' to 'right' over the
			plus strand amplicon) over the homology blocks and gaps. out[2] is True for gaps
			and False for homology blocks.
		"""
		c_s,c_q,c_r = self.block_sizes[0],self.q_starts[0],self.r_starts[0]
		### Output the first alignment block
		yield (int(c_q), int(c_q)+int(c_s)), (int(c_r), int(c_r)+int(c_s)), False
		for s,q,r in zip(self.block_sizes[1:],self.q_starts[1:],self.r_starts[1:]):
			### Output current gap
			yield (int(c_q)+int(c_s), int(q)), (int(c_r)+int(c_s), int(r)), True
			### Output current alignment block
			yield (int(q), int(q)+int(s)), (int(r), int(r)+int(s)), False
			c_s,c_q,c_r = s,q,r

	def get_alignment_blocks(self):
		"""
			Generator function that effectively iterates over all of the homology blocks of
			a given alignment.
		"""
		yield from (((q0,q1),(r0,r1)) for (q0,q1),(r0,r1),gap in self.get_alignment_groups() if not gap)

	def get_alignment_gaps(self):
		"""
			Generator function that iterates over all of the gaps of a given alignment.
		"""
		yield from (((q0,q1),(r0,r1)) for (q0,q1),(r0,r1),gap in self.get_alignment_groups() if gap)

	def calculate_score(self):
		gaps = sum([r1-r0 for (q0,q1),(r0,r1) in self.get_alignment_gaps()])
		return (self.matches - self.mismatches - gaps)/(self.q1-self.q0)

	def __str__(self):
		"""
			Lazy initialization of the string method.
		"""
		return self.string

def psl_iter(psl_path):
	"""
		Generator that iterates over all alignments in a psl file. 
		Returns a PslEntry object which has a rich API for dealing with 
		BLAT alignments.
	"""
	yield from (PslEntry(record) for record in tsv_iter(psl_path))

def psl_best(psl_path):
	"""
		Generator that iterates over the 'highest scoring' alignments for 
		each query sequence. 
		Returns a PslEntry object which has a rich API for dealing with 
		BLAT alignments.
	"""
	alignments = itertools.groupby(psl_iter(psl_path), lambda x: x.q_name)
	for q_name, alignment_group in alignments:
		#### Keeping only the best alignment
		for alignment in heapq.nlargest(1, alignment_group, key = lambda x: x.calculate_score()):
			yield alignment

def load_manifest(manifest_path):
	with open(manifest_path, "r") as manifest_file:
		manifest = yaml.load(manifest_file)
	return manifest

if __name__ == "__main__":
	psl_path = "/net/fileserv0/vol2/dchee7/genotyping/SMARCA/180226_MN00151_0152_A000H2GF3L/alignments/i7_1_N502/TL4407+TL4408_alignments.psl"
	for entry in itertools.islice(psl_iter(psl_path), 10):
		print(entry.count)
else:
	pass

