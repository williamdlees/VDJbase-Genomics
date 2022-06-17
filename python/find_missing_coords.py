# attempt to find missing genes in the franken by searching for the coding sequence
import os
from operator import itemgetter

from receptor_utils import simple_bio_seq as simple
import csv
from collections import namedtuple
import re

missing_genes = ['IGHV3-51-1', 'IGHV4-38-2', 'IGHV3-9', 'IGHV8-51-1', 'IGHV3-64D', 'IGHV3-25', 'IGHV7-77', 'IGHV1-38-4', 'IGHV1-69-2', 'IGHV3-NL1', 'IGHV4-30-1', 'IGHV3-62', 'IGHV1-8', 'IGHD5-24', 'IGHD1-1', 'IGHD6-25', 'IGHD1-20', 'IGHD6-19', 'IGHD5-12', 'IGHD5-5', 'IGHV3-64D', 'IGHD4-4', 'IGHD4-23', 'IGHD4-11', 'IGHD4-17', 'IGHD6-13', 'IGHD5-18', 'IGHD', 'IGHD6-6', 'IGHD2-21', 'IGHD1-7', 'IGHD1-14', 'IGHD1-26', 'IGHD7-27']
franken_file = '2021-12-08_franken_ref/igh.fasta'
imgt_files = [
    'gene_order/imgt_Homo sapiens_IGHV.fasta',
    'gene_order/imgt_Homo sapiens_IGHD.fasta',
    'gene_order/imgt_Homo sapiens_IGHJ.fasta',
]

# Directory containing the bed files to update
bed_dir = '2021-12-13_franken_gene_coords'

# L_Part1 sequences for V-genes we want. from IMGT000035 or deduced by hand using sequences from Mikocziova et al
imgt_lpart1 = {
    'IGHV4-38-2*02': 'atgaagcacctgtggtttttcctcctgctggtggcagctcccagat', # AC233755
    'IGHV3-9*01': 'atggagttgggactgagctggattttccttttggctattttaaaag', # M99651
    'IGHV8-51-1*02': 'atgcttgtatgtgtgcttttgtattctttcagattatttg',  # IMGT000035
    'IGHV3-64D*06': 'atggagttctggctgagctgggttctccttgttgccattttaaaag',  # IMGT000035
    'IGHV3-25*03': 'atggagtttgtgctgagctgggttttccttgttgctattttaaaac',  # IMGT000035 pseudogene? IMGT says stop codon in V-region
    'IGHV1-38-4*01': 'atggactggaactggaggatcctttttttggtggtggcaactacag', # KF698736
    'IGHV1-69-2*01': 'atggactgcacctggaggatcctcctcttggtggcagcagctacag', # IMGT000035
    'IGHV3-62*03': 'atggagtttgggctgagctgggttttccttgttgctattttaagag', # IMGT000035
    'IGHV1-8*01': 'atggactggacctggaggatcctcttcttggtggcagcagctacaa' # M99637,
}

imgt_lpart1 = {k: v.upper() for k, v in imgt_lpart1.items()}

# exclude_matches is used to omit duplicated annotations caused where the same sequence is found in the position for >1 gene

exclude_matches = {
    'IGHD5-5*01': 33042,
    'IGHD4-4*01': 43116,
    'IGHD4-11*01': 52722,
    'IGHD5-18*01': 51757
}

# Records to add to files

files_to_update = ['exon_1', 'exon_2', 'gencode_intron', 'genes', 'heptamer', 'nonamer', 'spacer']
added_records = {k: [] for k in files_to_update}


# Find the co-ordinates immediately below and above a position, given a sorted list of coords
def find_adjacent_coords(pos, coords):
    if pos <= coords[0]:
        return (None, coords[0])

    if pos >= coords[-1]:
        return(coords[-1], None)

    for i in range(len(coords) - 1):
        if coords[i] <= pos <= coords[i + 1]:
            return(coords[i], coords[i + 1])

    return (None, None)     # shouldn't happen


def print_v_recs(start, end, gene, allele):
    print(f"nonamer: igh {start-37} {start-29} {gene} nonamer {franken_3_to_5[start-37:start-29]}")
    added_records['nonamer'].append(['igh', start-37, start-29, gene, 'nonamer', franken_3_to_5[start-37:start-29]])

    print(f"spacer: igh {start-29} {start-7} {gene} spacer {franken_3_to_5[start-29:start-7]}")
    added_records['spacer'].append(['igh', start-29, start-7, gene, 'spacer', franken_3_to_5[start-29:start-7]])

    print(f"heptamer: igh {start-7} {start} {gene} heptamer {franken_3_to_5[start-7:start]}")
    added_records['heptamer'].append(['igh', start-7, start, gene, 'heptamer', franken_3_to_5[start-7:start]])

    print(f"exon_2: igh {start} {end+12} {gene}")
    added_records['exon_2'].append(['igh', start, end+12, gene])

    if allele in imgt_lpart1:
        target_lpart1 = simple.reverse_complement(imgt_lpart1[allele])
        for i in range(end+12+len(target_lpart1), end+12+len(target_lpart1) + 400):
            if franken_3_to_5[i:i+len(target_lpart1)] == target_lpart1:
                print(f"gencode_intron: igh {end+12} {i} {gene}")
                added_records['gencode_intron'].append(['igh', end+12, i, gene])

                print(f"l_part1: igh {i} {i + len(target_lpart1)} {gene}")
                added_records['exon_1'].append(['igh', i, i + len(target_lpart1), gene])

                added_records['genes'].append(['igh', start-37, end+12, gene])
    else:
        print('lpart1 not available')


def print_d_recs(start, end, gene):
    print(f"nonamer: igh {start-37} {start-29} {gene} nonamer {franken_3_to_5[start-37:start-29]}")
    added_records['nonamer'].append(['igh', start-37, start-29, gene, 'nonamer', franken_3_to_5[start-37:start-29]])

    print(f"spacer: igh {start-29} {start-7} {gene} spacer {franken_3_to_5[start-29:start-7]}")
    added_records['spacer'].append(['igh', start-29, start-7, gene, 'spacer', franken_3_to_5[start-29:start-7]])

    print(f"heptamer: igh {start-7} {start} {gene} heptamer {franken_3_to_5[start-7:start]}")
    added_records['heptamer'].append(['igh', start-7, start, gene, 'heptamer', franken_3_to_5[start-7:start]])

    print(f"exon_1: igh {start} {end} {gene}")
    added_records['exon_1'].append(['igh', start, end, gene])

    print(f"heptamer: igh {end} {end+7} {gene} heptamer {franken_3_to_5[end:end+7]}")
    added_records['heptamer'].append(['igh', end, end+7, gene, 'heptamer', franken_3_to_5[end:end+7]])

    print(f"spacer: igh {end+7} {end+30} {gene} spacer {franken_3_to_5[end+7:end+30]}")
    added_records['spacer'].append(['igh', end+7, end+30, gene, 'spacer', franken_3_to_5[end+7:end+30]])

    print(f"nonamer: igh {end+30} {end+38} {gene} nonamer {franken_3_to_5[end+30:end+38]}")
    added_records['nonamer'].append(['igh', end+30, end+38, gene, 'nonamer', franken_3_to_5[end+30:end+38]])

    added_records['genes'].append(['igh', start-37, end+38, gene])



franken_3_to_5 = simple.read_single_fasta(franken_file)
franken_5_to_3 = simple.reverse_complement(franken_3_to_5)
franken_len = len(franken_5_to_3)

imgt_seqs = {}
for imgt_file in imgt_files:
    seqs = simple.read_fasta(imgt_file)
    seqs = {k: v.replace('.', '') for k, v in seqs.items()}
    imgt_seqs.update(seqs)

gene_coords = {}
Coord = namedtuple('Coord', 'assembly start end gene')
with open('2021-12-13_franken_gene_coords/genes.bed', 'r') as fi:
    reader = csv.reader(fi, delimiter='\t')
    for row in reader:
        row = Coord(*row)
        gene_coords[int(row.start)] = row

sorted_starts = sorted(list(gene_coords.keys()))

for missing_gene in missing_genes:
    missing_alleles = [allele for allele in imgt_seqs.keys() if missing_gene == allele.split('*')[0]]
    if not missing_alleles:
        print(f"{missing_gene} has no alleles in the IMGT set")
        continue

    printed_header = False
    for missing_allele in missing_alleles:
        missing_seq = simple.reverse_complement(imgt_seqs[missing_allele])
        for match in re.finditer(missing_seq, franken_3_to_5):
            if missing_allele in exclude_matches and exclude_matches[missing_allele] == match.start():
                print(f"{missing_allele}: position {match.start()} excluded from consideration (duplicated gene sequence)")
                continue

            if not printed_header:
                print(f"allele {missing_allele} (5 to 3): {simple.reverse_complement(missing_seq)}\n")
                printed_header = True

            adjacent_below, adjacent_above = find_adjacent_coords(match.start(), sorted_starts)

            if adjacent_below:
                print(f"below: {gene_coords[adjacent_below].gene} ({gene_coords[adjacent_below].start}, {gene_coords[adjacent_below].end})")

            print(f"{missing_allele}: ({match.start()}, {match.end()}) (0-based, half-open) reverse: ({franken_len - match.end() + 1}, {franken_len - match.start()}) (1-based, closed)")

            if adjacent_above:
                print(f"above : {gene_coords[adjacent_above].gene} ({gene_coords[adjacent_above].start}, {gene_coords[adjacent_above].end})")

            if 'IGHV' in missing_gene:
                print_v_recs(match.start(), match.end(), missing_gene, missing_allele)
            elif 'IGHD' in missing_gene:
                print_d_recs(match.start(), match.end(), missing_gene)

        if not printed_header:
            print(f"{missing_allele}: no matches")
    print("\n\n")

    # Update files

for bedfile in files_to_update:
    oldfile = bedfile + '.bed.old'
    if os.path.isfile(oldfile):
        os.remove(oldfile)
    os.rename(os.path.join(bed_dir, bedfile + '.bed'), os.path.join(bed_dir, oldfile))

    with open(os.path.join(bed_dir, oldfile), 'r') as fi, open(os.path.join(bed_dir, bedfile + '.bed'), 'w', newline='') as fo:
        reader = csv.reader(fi, delimiter='\t')
        writer = csv.writer(fo, delimiter='\t')
        rows = [row for row in reader if row[1]]
        included_genes = [row[3] for row in rows]
        added = []
        for row in added_records[bedfile]:
            #if '64D' in row[3]:
            #    breakpoint()

            if row[3] in included_genes and row[1] not in added:
                print(f"replacing existing record for {row[3]} in {bedfile} ")
                for i in range(len(rows)):
                    if rows[len(rows) - i - 1][3] == row[3]:
                        del rows[len(rows) - i - 1]

            if row[1] not in added:
                rows.append(row)
                added.append(row[1])
            else:
                print(f"skipping duplicated new annotation for {row[3]}")
        for row in rows:
            row[1] = int(row[1])
            row[2] = int(row[2])
        rows.sort(key=itemgetter(1))
        writer.writerows(rows)
