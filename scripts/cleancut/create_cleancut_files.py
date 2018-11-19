#!/usr/bin/env python3

import csv
import re
import yaml

from collections import Counter

data_map = {
    "owner": 1,
    "desc": 2,
    "i7": 5,
    "i5": 6,
    "tl1": 8,
    "tl2": 9,
    "talen_sequence": 10,
    "talen_name": 12,
    "analysis": 13,
    "amplicon_sequence": 14,
    "hdr_amplicon_sequence": 15,
    "amplicon_name": 16,
}


def row_to_dict(row):
    d = {}
    for k, v in data_map.items():
        d[k] = row[v].strip()
        # N/A is a terrible sentinel value, replace with the empty string
        if d[k] == "N/A" :
            d[k] = ""

    # Remove spaces and uppercase the sequence
    for k in ["amplicon_sequence", "talen_sequence", "hdr_amplicon_sequence"]:
        d[k] = d[k].replace(" ", "").upper()

    # Pick various fallback names for people who won't name their talens (the fools)
    if not d["talen_name"]:
        if d['tl1'] and d['tl2']:
            d["talen_name"] = f"{d['tl1']}_{d['tl2']}"
        else:
            d["talen_name"] = f"{d['owner']}_talen"

    return d

def get_amplicons(data):
    amplicons = set((d['amplicon_name'], d['amplicon_sequence']) for d in data)
    return amplicons

def get_talens(data):
    talens = set((d['talen_name'], d['talen_sequence']) for d in data)
    return talens

def uniqueify_amplicons(data):
    amplicons = get_amplicons(data)
    # Sometimes people reuse names and we need to uniqueify them
    counts = Counter(a[0] for a in amplicons)
    for k, v in counts.items():
        if v == 1:
            continue
        to_update = set(a for a in amplicons if a[0] == k)
        amplicons = set(a for a in amplicons if a[0] != k)
        i = 1
        for a in to_update:
            new_name = f"{a[0]}_{i}"
            i += 1
            amplicons.add((new_name, a[1]))
            for j in range(len(data)):
                if data[j]['amplicon_sequence'] == a[1]:
                    data[j]['amplicon_name'] = new_name
    return data


def uniqueify_talens(data):
    talens = get_talens(data)
    # Sometimes people reuse names and we need to uniqueify them
    counts = Counter(t[0] for t in talens)
    for k, v in counts.items():
        if v == 1:
            continue
        to_update = set(t for t in talens if t[0] == k)
        talens = set(t for t in talens if t[0] != k)
        i = 1
        for t in to_update:
            new_name = f"{t[0]}_{i}"
            i += 1
            talens.add((new_name, t[1]))
            for j in range(len(data)):
                if data[j]['talen_sequence'] == t[1]:
                    data[j]['talen_name'] = new_name
    return data

def sanitize_data(data):
    r = re.compile(r"[_+/]+")
    for d in data:
        for k in "talen_name", "amplicon_name":
            if d[k] == "N/A":
                d[k] = "NA"
            d[k] = r.sub("_", d[k])
    return data


def read_data(filename):
    data = []
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        for row in reader:
            if row[0] == "FLOWCELL":
                continue
            d = row_to_dict(row)
            data.append(d)
    return data

def create_manifest(data):
    amplicons = get_amplicons(data)
    talens = get_talens(data)
    total_amplicons = amplicons | set([(d['amplicon_name'] + "_HDR", d['hdr_amplicon_sequence']) for d in data if d['hdr_amplicon_sequence']])
    return yaml.dump({
        "talens": {t[0]: t[1] for t in talens},
        "amplicons": {a[0]: a[1] for a in total_amplicons},
    }, default_flow_style=False)

def create_spec_sheet(data):
    return "\n".join(
        "\t".join([
            d['talen_name'],
            d['i7'],
            d['i5'],
            d['analysis'],
            d['amplicon_name'],
            d['amplicon_name'] + "_HDR" if d['hdr_amplicon_sequence'] else "NA"
        ]) for d in data)


def main(*args):
    data = read_data('source.tsv')
    data = sanitize_data(data)
    data = uniqueify_amplicons(data)
    data = uniqueify_talens(data)

    manifest = create_manifest(data)
    specsheet = create_spec_sheet(data)
    with open('spec-sheet.tsv', 'w') as f:
        f.write(specsheet)
    with open('manifest.yaml', 'w') as f:
        f.write(manifest)


if __name__ == "__main__":
    main()

