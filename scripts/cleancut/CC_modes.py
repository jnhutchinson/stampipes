### Copyright Daniel R. Chee 2017
### Altius Institute of Biomedical Sciences

import argparse
import collections
import itertools
import math
import numpy
import pyfaidx
import re
import regex

import microhomology
import TalenModel
import file_parsers

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#def dimer_efficiency(dimer_id, talen_path, read_fa, amplicon_fa, alignment_path, prefix="./test"):
def dimer_efficiency(dimer_id, talen_path, read_fa, amplicon_fa, alignment_path, prefix="./test"):
    ###
    talens = file_parsers.load_manifest(talen_path)["talens"]
    ###
    target_str = ""
    if dimer_id in talens:
        target_str = f"{talens[dimer_id]}"
    else:
        for monomer in dimer_id.split("+"):
            if monomer in talens:
                target_str = f"{target_str}[ATCGN]*{talens[monomer]}" if target_str else f"{talens[monomer]}"
    #talen_re = re.compile(f"{talens[left_monomer]}[ATCGN]*{talens[right_monomer]}")
    talen_re = re.compile(f"{target_str}") if target_str else re.compile(f"[ATCGN]*")
    ###
    cutsite = {}
    amplicon_sequence = {}
    amplicons = pyfaidx.Fasta(amplicon_fa)
    for key in amplicons.keys():
        sequence = str(amplicons[key])
        cutsite[key] = talen_re.search(sequence).span() if talen_re.search(sequence) else (0,len(sequence))
        amplicon_sequence[key] = sequence
        if key == "WT":
            length = len(sequence)
    ###
    indels = collections.Counter()
    mismatches = numpy.zeros(length)
    dimer_total = 0
    dimer_seen = set()
    total_reads = 0
    ### 
    read_fasta = pyfaidx.Fasta(read_fa)
    ### First pass builds the mutational landscape
    regex = re.compile(r"read_[0-9]*_[0-9]*")
    for record in itertools.islice((x for x in file_parsers.psl_iter(alignment_path) if x.r_name == "WT"), None):
        ### Number of reads that this read corresponds to. 
        for (q0,q1),(r0,r1),gap in record.get_alignment_groups():
            if not gap:
                ### Check for mismatches
                for r_i,ref_base,query_base in zip(range(r0,r1), amplicon_sequence[record.r_name][r0:r1], \
                        str(read_fasta[record.q_name][q0:q1])):
                    ###
                    if ref_base != query_base:
                        mismatches[r_i]+=record.count
                #coverage[r0:r1]+=record.count
            else:
                ### Aligning the read
                ### Checking for overlap with the .....
                if (r0 <= cutsite[record.r_name][1]) and (cutsite[record.r_name][0] <= r1):
                    indels[(r0,r1)]+=record.count
                    if record.q_name not in dimer_seen:
                        dimer_total+=record.count
                        dimer_seen.add(record.q_name)
                ###
        total_reads+=record.count

    snp_thresh = numpy.mean(mismatches/total_reads) + (2*(numpy.mean(mismatches/total_reads))) if total_reads else 0
    #snp_thresh = 0.01

    ### Second pass outputs the alignments
    HDR_total = 0
    aligned_seqs = collections.OrderedDict()
    for record in itertools.islice(file_parsers.psl_best(alignment_path), None):
        offset = 0
        diff = 0
        sequence = []
        for (q0,q1),(r0,r1),gap in record.get_alignment_groups():
            ### adding the beginning sequence if needed
            if len(sequence) == 0 and r0 > 0:
                sequence = sequence + [base for base in amplicon_sequence[record.r_name][:r0]]
            if not gap:
                ### Check for mismatches
                for r_i,ref_base,query_base in zip(range(r0,r1), amplicon_sequence[record.r_name][r0:r1], \
                        str(read_fasta[record.q_name][q0:q1])):
                    ###
                    if ref_base != query_base:
                        if record.r_name == "WT":
                            #if mismatches[r_i]/total_reads >= 0.01:
                            if mismatches[r_i]/total_reads >= snp_thresh:
                                ### Actual mismatch
                                sequence.append(query_base.lower())
                            else:
                                sequence.append(ref_base)
                        else:
                            sequence.append(query_base.lower())
                    else:
                        sequence.append(query_base)
            else:
                if q1-q0 > 0:
                    ### Small insertion
                    sequence.append(str(read_fasta[record.q_name][q0:q1]).lower())
                    diff+=q1-q0
                ### Aligning the read
                if r1-r0 > 0:
                    sequence.append("-"*(r1-r0))
                    diff-=r1-r0
        if record.r_name == "donor":
            HDR_total += record.count
        ### Assume the end of the reads is "wildtype"
        sequence = sequence + [base for base in amplicon_sequence[record.r_name][r1:]]
        ###
        sequence = "".join(sequence)
        if (record.r_name,diff,sequence) not in aligned_seqs:
            aligned_seqs[(record.r_name,diff,sequence)] = 0
        aligned_seqs[(record.r_name,diff,sequence)] += record.count
    ### Creating a list of indel size ranges
    indel_sizes = (list(range(1,11)) + [15, 20, math.inf])[::-1]
    indel_dist = numpy.zeros(length)
    deletions = {}
    cleavages = {}
    ### Creating a dictionary of numpy arrays.
    for i in indel_sizes:
        deletions[f"<={i}"] = numpy.zeros(length)
        cleavages[f"<={i}"] = numpy.zeros(length)
    ### Getting a set of predicted microhomology
    WT_seq = amplicon_sequence["WT"]
    micro_total = 0
    ### Iterating over indels and ... 
    for (r0,r1),count in indels.items():
        indel_dist[r1-r0]+=count
        for i in indel_sizes:
            if (r1-r0) <= i:
                deletions[f"<={i}"][r0:r1]+=count
                cleavages[f"<={i}"][r0]+=count
                cleavages[f"<={i}"][r1]+=count
            else:
                break
        m = 0
        while ((r0+m) < r1) and (WT_seq[r0+m] == WT_seq[r1+m]):
            m+=1
        micro_total += count if m >= 3 else 0
    ### Outputting the alignments.
    with open(f"{prefix}_aligned.txt", "w") as aligned_file:
        for (geno,diff,seq),count in sorted(aligned_seqs.items(), key=lambda x: x[1], reverse=True):
            print(geno,diff,seq,count,f"{(count/total_reads)*100:.2f}%", sep="\t",file=aligned_file)
    ### Writing the efficiency file.
    effic_path = f"{prefix}_efficiency.txt"
    with open(effic_path, "w") as effic_file:
        print("dimer(s)\tcut\tHDR\tmicro\ttotal\tcut %\tHDR %\tmicro %", file=effic_file)
        dimer_effic = round((dimer_total/total_reads)*100, 2) if total_reads else 0
        HDR_effic = round((HDR_total/total_reads)*100, 2) if total_reads else 0
        micro_effic = round((micro_total/total_reads)*100, 2) if total_reads else 0
        print(f"{dimer_id}\t{dimer_total}\t{HDR_total}\t{micro_total}\t{total_reads}\t{dimer_effic}\t{HDR_effic}\t{micro_effic}", file=effic_file)
    ### Outputting the deletion and cleavage profiles
    deletion_path = f"{prefix}_deletion_profile.txt"
    cleavage_path = f"{prefix}_cleavage_profile.txt"
    with open(deletion_path, "w") as deletion_file, open(cleavage_path, "w") as cleavage_file:
        ###
        cols = ["position", "base"]
        for i in indel_sizes[::-1]:
            cols.append(f"<={i}")
            cols.append(f"<={i}%")
        print(*cols, sep="\t",file=deletion_file)
        print(*cols, sep="\t",file=cleavage_file)
        ###
        for position in range(length):
            delete = [position, amplicon_sequence["WT"][position]] 
            cleave = [position, amplicon_sequence["WT"][position]] 
            for i in indel_sizes[::-1]:
                delete.append(deletions[f"<={i}"][position])
                delete.append(deletions[f"<={i}"][position]/total_reads)
                cleave.append(cleavages[f"<={i}"][position])
                cleave.append(cleavages[f"<={i}"][position]/total_reads)
            print(*delete, sep="\t",file=deletion_file)
            print(*cleave, sep="\t",file=cleavage_file)
    ### Storing the indel distribution
    numpy.save(f"{prefix}_indel_dist.npy", indel_dist)
    ### Outputting the stacked deletion plots
    plt.figure(figsize=(10,8))
    plt.suptitle(f"{dimer_id} Stacks")
    gs = gridspec.GridSpec(len(indel_sizes), 1, hspace=0.35)
    prev = None
    for n,i in enumerate(indel_sizes[::-1]):
        if prev:
            prev.set_xticks([])
            prev.set_xticklabels([])
        ax = plt.subplot(gs[n, 0])
        ax.plot(deletions[f"<={i}"]/total_reads)
        ax.fill_between(numpy.arange(0, len(deletions[f"<={i}"])), 0, deletions[f"<={i}"]/total_reads, alpha=0.8)
        max_val = numpy.max(deletions[f"<={math.inf}"])/total_reads
        ### Axis on the right side
        ax1 = ax.twinx()
        ax1.set_yticks([0, .9])
        ax1.set_yticklabels(["0.000", f"{max_val:.2}"])
        for tick in ax1.yaxis.get_major_ticks():
            tick.label.set_fontsize(4) 
        ###
        ax.set_ylim(0,max_val + (max_val*.1))
        ax.set_yticks([(max_val + (max_val*.1))/2])
        ax.set_yticklabels([f"<={i}"])
        prev = ax
    plt.savefig(f"{prefix}_stacked_deletion.pdf")
    ###
    def plot_talens(ax):
        ax.axvspan(cutsite["WT"][0],cutsite["WT"][1], color="#FFBB00", alpha=0.4)
        return
    ###
    plt.figure(figsize=(10,8))
    gs = gridspec.GridSpec(3, 1)
    gs.update(hspace=0.5)
    ### INDEL information
    ax = plt.subplot(gs[0, 0])
    ax.set_title("FOK1 Cleavages")
    ax.plot(cleavages["<=inf"]/total_reads)
    ax.set_ylim(0,1.05)
    plot_talens(ax)
    ## Deletion Plot
    ax = plt.subplot(gs[1, 0])
    ax.set_title("Deleted Bases")
    ax.plot(deletions["<=inf"]/total_reads)
    ax.set_ylim(0,1.05)
    plot_talens(ax)
    ### Indel Size distribution
    ax = plt.subplot(gs[2, 0])
    ax.set_title("Indel Size Distribution")
    ax.bar(numpy.arange(len(indel_dist)), indel_dist/numpy.max(indel_dist))
    ax.set_ylim(0,1.05)
    ###
    plt.savefig(f"{prefix}_dimer_plot.pdf")
    return

def cinco_de_mayo(dimer_id, talen_path, read_fa, amplicon_fa, alignment_path, prefix="./test"):

    ###
    talens = file_parsers.load_manifest(talen_path)["talens"]
    ###
    target_str = ""
    if dimer_id in talens:
        target_str = f"{talens[dimer_id]}"
    else:
        for monomer in dimer_id.split("+"):
            if monomer in talens:
                target_str = f"{target_str}[ATCGN]*{talens[monomer]}" if target_str else f"{talens[monomer]}"

    talen_re = re.compile(f"{target_str}") if target_str else re.compile(f"[ATCGN]*")

    cutsite = {}
    amplicon_sequence = {}
    amplicons = pyfaidx.Fasta(amplicon_fa)
    for key in amplicons.keys():
        sequence = str(amplicons[key])
        cutsite[key] = talen_re.search(sequence).span() if talen_re.search(sequence) else (0,len(sequence))
        amplicon_sequence[key] = sequence
        if key == "WT":
            length = len(sequence)

    indels = collections.Counter()
    mismatches = collections.defaultdict(int)
    dimer_total = 0
    dimer_seen = set()
    total_reads = 0

    regex = re.compile(r"read_[0-9]*_[0-9]*")

    # Building the sequencing error landscape
    read_fasta = pyfaidx.Fasta(read_fa)
    for record in file_parsers.psl_iter(alignment_path):

        if record.r_name == "WT":
            # Number of reads to which this alignment corresponds.
            if regex.match(record.q_name):
                inc = int(record.q_name.split("_")[-1]) 
            else:
                inc = 1

            for (q0,q1),(r0,r1),gap in record.get_alignment_groups():
                if not gap:
                    ### Check for mismatches

                    align_zip = zip(
                        range(r0,r1),
                        amplicon_sequence[record.r_name][r0:r1],
                        str(read_fasta[record.q_name][q0:q1])
                    )

                    for r_i,ref_base,query_base in align_zip:
                        if ref_base != query_base:
                            mismatches[(r_i,query_base)]+=inc

            total_reads+=inc

    # Second pass counts alleles and aligns the sequences.
    allele_count = collections.defaultdict(int)
    allele_sequence = {}
    for record in file_parsers.psl_best(alignment_path):

        #
        mutational_profile = [record.r_name]

        # Number of reads that this read corresponds to. 
        if regex.match(record.q_name):
            inc = int(record.q_name.split("_")[-1])
        else:
            inc = 1

        diff = 0
        first = True
        sequence = []
        for (q0,q1),(r0,r1),gap in record.get_alignment_groups():

            ### adding the beginning sequence if needed
            if len(sequence) == 0 and r0 > 0:
                sequence = sequence + [base for base in amplicon_sequence[record.r_name][:r0]]

            if not gap:
                ### Check for mismatches
                for r_i,ref_base,query_base in zip(range(r0,r1), amplicon_sequence[record.r_name][r0:r1], str(read_fasta[record.q_name][q0:q1])):
                    if ref_base != query_base:
                        if mismatches[r_i]/total_reads >= 0.01:
                            ### Actual mismatch
                            mutational_profile.append(
                                (
                                    "SNP",
                                    r_i,
                                    query_base
                                )
                            )
                            sequence.append(query_base.lower())
                        else:
                            sequence.append(ref_base)
                    else:
                        sequence.append(query_base)
            else:
                ### Aligning the read
                if q1-q0 > 0:

                    # insertion
                    insertion = str(read_fasta[record.q_name][q0:q1]).lower()

                    mutational_profile.append(
                        (
                            "insertion",
                            r0,
                            r1,
                            insertion
                        )
                    )
                    sequence.append(insertion)
                    diff+=q1-q0
                if r1-r0 > 0:
                    mutational_profile.append(
                        (
                            "deletion",
                            r0,
                            r1
                        )
                    )
                    sequence.append("-"*(r1-r0))
                    diff-=r1-r0
        
        if r1 < len(amplicon_sequence[record.r_name])-1:
            sequence.append(amplicon_sequence[record.r_name][r1:])

        allele_count[tuple(mutational_profile)]+=inc
        allele_sequence[tuple(mutational_profile)] = "".join(sequence)

    ### 
    total = 0
    prominant_alleles = []
    for x,c in allele_count.items():
        if c/total_reads >= 0.01:
            prominant_alleles.append((x,c))
            total+=c

    geno_path = f"{prefix}_genotypes.txt"
    with open(geno_path, "w") as geno_file:
        print(
            "dimer",
            "genotype",
            "call",
            "sequence",
            "count",
            "read %",
            sep="\t",
            file=geno_file
        )

        ### First pass builds the genotype string
        genotype = set()
        for profile,count in prominant_alleles:
            
            call = profile[0]            

            if len(profile) > 1:
                genotype.add(f"{call}*")
            else:
                genotype.add(f"{call}")

        genotype="/".join(sorted(genotype, reverse=True))

        ### Second pass prints the output
        for profile,count in prominant_alleles:

            call = profile[0]            

            if len(profile) > 1:
                call = f"{call}*"

            print(
                dimer_id,
                genotype,
                call,
                allele_sequence[profile],
                count,
                f"{(count/total_reads)*100:.0f}%",
                sep="\t",
                file=geno_file
            )

    return

def setup_parser():

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    
    dimer_efficiency_parser = subparsers.add_parser("dimer-efficiency")
    dimer_efficiency_parser.add_argument("dimer_id")
    dimer_efficiency_parser.add_argument("talen_path")
    dimer_efficiency_parser.add_argument("read_fa")
    dimer_efficiency_parser.add_argument("amplicon_fa")
    dimer_efficiency_parser.add_argument("alignment_path")
    dimer_efficiency_parser.add_argument("--prefix", type=str, default="./test")

    cinco_de_mayo_parser = subparsers.add_parser("cinco-de-mayo")
    cinco_de_mayo_parser.add_argument("dimer_id")
    cinco_de_mayo_parser.add_argument("talen_path")
    cinco_de_mayo_parser.add_argument("read_fa")
    cinco_de_mayo_parser.add_argument("amplicon_fa")
    cinco_de_mayo_parser.add_argument("alignment_path")
    cinco_de_mayo_parser.add_argument("--prefix", type=str, default="./test")

    return parser.parse_args()

def run(args):

    if args.command == "dimer-efficiency":

        dimer_efficiency(
            args.dimer_id,
            args.talen_path, 
            args.read_fa, 
            args.amplicon_fa, 
            args.alignment_path,
            prefix=args.prefix
        )

    if args.command == "cinco-de-mayo":

        cinco_de_mayo(
            args.dimer_id,
            args.talen_path, 
            args.read_fa, 
            args.amplicon_fa, 
            args.alignment_path, 
            prefix=args.prefix
        )

    return

if __name__ == "__main__":
    run(setup_parser())

