# make GENE file

import os
import csv
from operator import attrgetter, itemgetter

coords = {}

with open('NONAMER.bed', 'r') as fi:
    reader = csv.DictReader(fi, delimiter='\t', fieldnames=['chain', 'start', 'end', 'gene'])
    for row in reader:
        rec = {
            'chain': 'igl',
            'start': 0,
            'end': int(row['end']),
            'gene': row['gene'],
        }
        coords[row['gene']] = rec


with open('UTR.bed', 'r') as fi:
    reader = csv.DictReader(fi, delimiter='\t', fieldnames=['chain', 'start', 'end', 'gene'])
    for row in reader:
        coords[row['gene']]['start'] = int(row['start'])

recs = list(coords.values())
recs = sorted(recs, key=itemgetter('start'))

with open('GENE.bed', 'w') as fo:
    for rec in recs:
        fo.write(f"{rec['chain']}\t{rec['start']}\t{rec['end']}\t{rec['gene']}\n")

# EXON_2

coords = {}

with open('L-PART2.bed', 'r') as fi:
    reader = csv.DictReader(fi, delimiter='\t', fieldnames=['chain', 'start', 'end', 'gene'])
    for row in reader:
        rec = {
            'chain': 'igl',
            'start': int(row['start']),
            'end': 0,
            'gene': row['gene'],
        }
        coords[row['gene']] = rec

with open('HEPTAMER.bed', 'r') as fi:
    reader = csv.DictReader(fi, delimiter='\t', fieldnames=['chain', 'start', 'end', 'gene'])
    for row in reader:
        if 'V' in row['gene']:
            coords[row['gene']]['end'] = int(row['start'])

recs = list(coords.values())
recs = sorted(recs, key=itemgetter('start'))

with open('EXON_2.bed', 'w') as fo:
    for rec in recs:
        fo.write(f"{rec['chain']}\t{rec['start']}\t{rec['end']}\t{rec['gene']}\n")
