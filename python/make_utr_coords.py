# Make UTR co-ordinates using upstream sequences from Huang, Thornqvist and Ohlin

from receptor_utils import simple_bio_seq as simple
import csv
from operator import itemgetter

def read_bed(bedfile):
    with open(bedfile, 'r') as fi:
        genes = {}
        reader = csv.DictReader(fi, delimiter='\t', fieldnames=['chain', 'start', 'end', 'gene'])
        for row in reader:
            row['start'] = int(row['start'])
            row['end'] = int(row['end'])
            genes[row['gene']] = row

        return genes


exon_1 = read_bed('EXON_1.bed')
intron = read_bed('INTRON.bed')
upstream = simple.read_fasta('ref/Huang_Thornqvist_Ohlin_Supplementary_Data_2.txt')

processed_genes = []
utrs = []

for allele in upstream.keys():
    gene = allele.split('*')[0]

    if gene in processed_genes:
        continue

    processed_genes.append(gene)

    if gene not in exon_1 or gene not in intron:
        print(f'gene {gene} not found in bed files')
        continue

    if gene == 'IGHV4-34':
        breakpoint()

    leader_length = exon_1[gene]['end'] - exon_1[gene]['start'] + 11
    leader_seq = upstream[allele][-leader_length:]
    utr_seq = upstream[allele][:-leader_length]

    if leader_seq[:3] != 'ATG':
        print(f'Non-canonical leader for gene {gene}')

    utrs.append([
        'igh',
        str(exon_1[gene]['end']),
        str(exon_1[gene]['end'] + len(utr_seq)),
        gene,
    ])

unprocessed_genes = set(exon_1.keys()) - set(processed_genes)
unprocessed_genes = [x for x in unprocessed_genes if 'IGHV' in x]

if len(unprocessed_genes):
    print(f'Unprocessed: {unprocessed_genes}')

with open('UTR.bed_mod', 'w') as fo:
    for utr in sorted(utrs, key=lambda x: int(x[2])):
        fo.write('\t'.join(utr) + '\n')

