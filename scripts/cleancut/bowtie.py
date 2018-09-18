### Copyright Daniel R. Chee 2017
### Altius Institute of Biomedical Sciences

import csv
import tempfile
import itertools
import subprocess
import collections
import os
import sys
import json

config_file = os.environ.get("TALEN_CONFIG", None)

def get_tmpdir(config_file):
	with open(config_file) as raw_config:
		config = json.loads(raw_config.read())
		try:
			return config["tmpdir"]
		except KeyError:
			return "/tmp"
	raise Exception ("Cannot load config file {}".format(config_file))

def load_indices(config_file):
	with open(config_file) as raw_config:
		config = json.loads(raw_config.read())
		try:
			return config["bowtie_index"]
		except IndexError:
			raise Exception("No bowtie indexes specified (bowtie_index in config)")
	raise Exception ("Cannot load config file {}".format(config_file))

indices = load_indices(config_file)

def get_number_of_hits(sequence_dict, mismatches=0,assembly="hg38"):
	"""
		Function that takes a dictionary of id,sequence pairs
		and returns the number of genomic hits for each id.
	"""
	### Pointing to the location of my index genomes (indexed via bowtie command)
	try:
		 index = indices[assembly]
	except Error as e:
		print(f"Cannot find assembly {assembly}")
		sys.exit(1)
	### Letting the context manager take care of opening and closing the tempfiles
	with tempfile.NamedTemporaryFile(mode="w+t", dir=get_tmpdir(config_file), suffix=".fa") as fasta, tempfile.TemporaryFile(mode="w+t", dir=get_tmpdir(config_file), suffix=".bowtie") as alignments:
		### Loading sequence_dict into fasta format
		for seq_id,sequence in sequence_dict.items():
			print(f">{seq_id}\n{sequence}", file=fasta)
		fasta.seek(0)
		### Calling bowtie
		bowtie_c = ["bowtie","--quiet","--all","-v",f"{mismatches}",f"{index}","-f",f"{fasta.name}"]
		bowtie = subprocess.Popen(bowtie_c, stdout=alignments)
		bowtie.communicate()
		alignments.seek(0)
		### Parsing the number of hits simply by counting the number of times 
		### the seq_id appears in the output file.
		hits = (hit[0] for hit in csv.reader(alignments, delimiter='\t'))
		num_hits = collections.Counter()
		num_hits.update(hits)

	return num_hits

def get_coordinates(sequence_dict, mismatches=0,assembly="hg38"):
	### Pointing to the location of my index genomes (indexed via bowtie command)
	index = f"/home/audrakj/datastore/bowtie/{assembly}"
	### Letting the context manager take care of opening and closing the tempfiles
	with tempfile.NamedTemporaryFile(mode="w+t") as fasta, tempfile.TemporaryFile(mode="w+t") as alignments:
		### Loading sequence_dict into fasta format
		for seq_id,sequence in sequence_dict.items():
			print(f">{seq_id}\n{sequence}", file=fasta)
		fasta.seek(0)
		### Calling bowtie
		bowtie_c = ["bowtie","--quiet","--all","-v",f"{mismatches}",f"{index}","-f",f"{fasta.name}"]
		bowtie = subprocess.Popen(bowtie_c, stdout=alignments)
		bowtie.communicate()
		alignments.seek(0)
		### Parses through the alignments and returns the first genomic match. 
		### Keep this in mind if there are multiple genomic matches.
		matches = {}
		for line in csv.reader(alignments, delimiter='\t'):
			if line[0] not in matches:
				matches[line[0]] = (line[2], int(line[3]), int(line[3])+len(line[4]), line[1])
	return matches
