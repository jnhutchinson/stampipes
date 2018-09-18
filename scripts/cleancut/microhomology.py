### Copyright Daniel R. Chee 2018
### Altius Institute of Biomedical Sciences

import argparse
import regex

def mmej_iter(sequence, m_min=5, m_max=25):
	micro_re = regex.compile(f"(?P<micro>[ATCG]{{{m_min},{m_max}}})[ATCGN]*(?P=micro)")
	for match in micro_re.finditer(sequence,overlapped=True):
		x,y = match.span()
		for i in range(len(match.group("micro"))):
			yield x+i,y-(len(match.group("micro"))-i)
		#yield x,y-len(match.group("micro"))
		#yield x+len(match.group("micro")),y

def mmej_set(sequence, m_min=5, m_max=25):
	micro_set = set()
	for mmej in mmej_iter(sequence, m_min=m_min, m_max=m_max):
		micro_set.add(mmej)
	return micro_set

def display_microhomology(sequence, m_min=5, m_max=25):
	### Printing the wild type sequence
	print(sequence)
	micro_re = regex.compile(f"(?P<micro>[ATCG]{{{m_min},{m_max}}})[ATCGN]*(?P=micro)")
	for match in micro_re.finditer(sequence,overlapped=True):
		micro = []
		m0,m1 = match.span()
		micro.append(" "*(m0-0))
		micro.append(match.group("micro"))
		micro.append("."*((m1-len(match.group("micro")))-(m0+len(match.group("micro")))))
		micro.append(match.group("micro"))
		micro.append(" "*(len(sequence)-m1))
		print("".join(micro))

if __name__ == "__main__":
	###
	parser = argparse.ArgumentParser()
	###
	parser.add_argument("sequence")
	parser.add_argument("--m-min", type=int, default=5)
	parser.add_argument("--m-max", type=int, default=25)
	###
	args = parser.parse_args()
	###
	display_microhomology(args.sequence, m_min=args.m_min, m_max=args.m_max)
else:
	pass
