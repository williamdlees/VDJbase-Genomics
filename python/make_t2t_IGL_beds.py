# Make bed files for T2T IGL assembly using existing beds as a template

from receptor_utils import simple_bio_seq as simple
import csv

digger_file = '../t2t_coords/CHM13Y_IGL.csv'
old_utr_bed_file = 'pre_t2t/UTR.bed'
assembly_file = '../t2t_coords/CHM13Y_IGL.fa'
imgt_ref_file = '../t2t_coords/Homo_sapiens_IGLVJ.fasta'
locus = 'IGL'

# Offset to add to co-ordinates in digger file, eg to map to entire chromosome assembly
coord_offset = 0

beds = {
    'EXON_1': [],
    'EXON_2': [],
    'GENE': [],
    'HEPTAMER': [],
    'INTRON': [],
    'L-PART2': [],
    'NONAMER': [],
    'SPACER': [],
    'UTR': [],
}

# Read digger file
digger_recs = simple.read_csv(digger_file)

# read assembly file
assembly = simple.read_fasta(assembly_file)
assembly_name = list(assembly.keys())[0]
assembly = list(assembly.values())[0]

# Read old UTR bed file
utr_lengths = {}

# Read imgt reference
imgt_ref = simple.read_fasta(imgt_ref_file)

with open(old_utr_bed_file, 'r') as fi:
    reader = csv.DictReader(fi, delimiter='\t', fieldnames=['assembly', 'start', 'end', 'gene'])
    for row in reader:
        utr_lengths[row['gene'].replace('_5-UTR', '')] = int(row['end']) - int(row['start'])


# convert 1-based co-ords to 0-based, converting to int and taking account of coord_offset
def convert_coords(start, end):
    return int(start) - 1 + coord_offset, int(end) + coord_offset


# Return a sequence defined by bed co_ordinates
def get_seq(bed_name, gene):
    for rec in beds[bed_name]:
        if rec[2] == gene:
            start, end = rec[:2]
            return assembly[start:end]

    return None

# process functional records in digger file
for rec in digger_recs:
    # start positions 138274 (4-60), 732949 (3-16) are non-functional in T2T but known to be functional in other samples
    # start position 64630 (X-XX) 319484 (9-49) 579832 (3-32) leader has been mis-identified by digger
    if rec['functional'] == 'Functional' or rec['start'] in ('138274', '579832', '732949'):
        if rec['start'] == '579832':
            rec['start'] = int(rec['start']) + 20
            rec['l_part2_end'] = int(rec['l_part2_end']) + 20
            rec['l_part2_start'] = rec['l_part2_end'] - 10
            rec['l_part1_end'] = rec['l_part2_start'] - 147
            rec['l_part1_start'] = rec['l_part1_end'] - 42

        # as-yet unnamed apparently functional gene at this position
        if rec['start'] == '64630':
            rec['imgt_match'] = 'IGLVX-XX*01'

        if rec['start'] == '319484' or rec['start'] == '64630':
            rec['start'] = int(rec['start']) + 6
            rec['l_part2_end'] = int(rec['l_part2_end']) + 6

        if rec['imgt_match'] in imgt_ref:
            rec['imgt_seq'] = imgt_ref[rec['imgt_match']]
        else:
            rec['imgt_seq'] = ''

        gene = rec['imgt_match'].split('*')[0]

        if 'V' in rec['gene_type']:
            beds['EXON_1'].append((*convert_coords(rec['l_part1_start'], rec['l_part1_end']), gene))
            beds['EXON_2'].append((*convert_coords(rec['l_part2_start'], rec['end']), gene))
            beds['L-PART2'].append((*convert_coords(rec['l_part2_start'], rec['l_part2_end']), gene))
            beds['HEPTAMER'].append((*convert_coords(rec['3_rss_start'], int(rec['3_rss_start']) + 6), gene))
            beds['SPACER'].append((*convert_coords(int(rec['3_rss_start']) + 7, int(rec['3_rss_end']) - 9), gene))
            beds['NONAMER'].append((*convert_coords(int(rec['3_rss_end']) - 8, rec['3_rss_end']), gene))
            beds['INTRON'].append((*convert_coords(int(rec['l_part1_end']) + 1, int(rec['l_part2_start']) - 1), gene))

            if gene in utr_lengths and utr_lengths[gene] > 0:
                beds['UTR'].append((*convert_coords(int(rec['l_part1_start']) - utr_lengths[gene], int(rec['l_part1_start']) - 1), gene))
                beds['GENE'].append((*convert_coords(int(rec['l_part1_start']) - utr_lengths[gene], int(rec['3_rss_end'])), gene))
            else:
                beds['GENE'].append((*convert_coords(int(rec['l_part1_start']), int(rec['3_rss_end'])), gene))


            if(get_seq('EXON_1', gene) != rec['l_part1']):
                print(f"EXON_1 mismatch for {gene}")
            if(get_seq('L-PART2', gene) != rec['l_part2']):
                print(f"L-PART2 mismatch for {gene}")
            if(get_seq('HEPTAMER', gene) != rec['v_heptamer']):
                print(f"HEPTAMER mismatch for {gene}")
            if(get_seq('NONAMER', gene) != rec['v_nonamer']):
                print(f"NONAMER mismatch for {gene}")
            if(get_seq('EXON_2', gene)[11:] != rec['seq']):
                print(f"V-REGION mismatch for {gene}")
            if(rec['imgt_seq'] and get_seq('EXON_2', gene)[11:] != rec['imgt_seq']):
                print(f"Sequence mismatch with IMGT reference for {gene}")

        elif 'J' in rec['gene_type']:
            beds['HEPTAMER'].append((*convert_coords(int(rec['5_rss_end']) - 6, rec['5_rss_end']), gene))
            beds['SPACER'].append((*convert_coords(int(rec['5_rss_start']) + 9, int(rec['5_rss_end']) - 7), gene))
            beds['NONAMER'].append((*convert_coords(rec['5_rss_start'], int(rec['5_rss_start']) + 8), gene))
            beds['GENE'].append((*convert_coords(int(rec['5_rss_start']), int(rec['end'])), gene))
            beds['EXON_1'].append((*convert_coords(int(rec['start']), int(rec['end'])), gene))

            if(get_seq('HEPTAMER', gene) != rec['j_heptamer']):
                print(f"HEPTAMER mismatch for {gene}")
            if(get_seq('NONAMER', gene) != rec['j_nonamer']):
                print(f"NONAMER mismatch for {gene}")
            if(get_seq('EXON_1', gene) != rec['seq']):
                print(f"EXON_1 mismatch for {gene}")
            if(rec['imgt_seq'] and get_seq('EXON_1', gene) != rec['imgt_seq']):
                print(f"Sequence mismatch with IMGT reference for {gene}")


'''
The output checks the V-REGION or J-REGION against the IMGT reference.
In addition it checks most fields in the bed files against the sequences in Digger output
The fixes made at the start of record processing will cause mismatches relative to Digger output.

Expected output:
L-PART2 mismatch for IGLVX-XX
V-REGION mismatch for IGLVX-XX ...L_PART2 was not determined correctly bu Digger
Sequence mismatch with IMGT reference for IGLV4-60 ...T2T has a deletion in this gene
L-PART2 mismatch for IGLV9-49
V-REGION mismatch for IGLV9-49 ...L_PART2 was not determined correctly bu Digger
EXON_1 mismatch for IGLV3-32
L-PART2 mismatch for IGLV3-32
V-REGION mismatch for IGLV3-32 ... L_PART1 was not correctly determined by Digger
Sequence mismatch with IMGT reference for IGLV2-14  T2T has a novel allele for this gene: IGLV2-14*01_g168t
'''

# write bed files
for bed_name in beds:
    with open(f'{bed_name}.bed', 'w', newline='\n') as fo:
        writer = csv.writer(fo, delimiter='\t')
        writer.writerows([list((assembly_name, *row))for row in beds[bed_name]])

