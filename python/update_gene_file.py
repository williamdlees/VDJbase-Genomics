# update v-gene co-ordinates to include UTR

import csv


def read_bed(bedfile):
    with open(bedfile, 'r') as fi:
        genes = {}
        reader = csv.DictReader(fi, delimiter='\t', fieldnames=['chain', 'start', 'end', 'gene'])
        for row in reader:
            row['start'] = int(row['start'])
            row['end'] = int(row['end'])
            genes[row['gene']] = row

        return genes


genes = read_bed('GENE.bed')
utrs = read_bed('UTR.bed')

for gene, rec in utrs.items():
    if gene in genes:
        genes[gene]['end'] = utrs[gene]['end']

with open('GENE_new.bed', 'w') as fo:
    for gene, rec in genes.items():
        fo.write(f"{rec['chain']}\t{rec['start']}\t{rec['end']}\t{rec['gene']}\n")

