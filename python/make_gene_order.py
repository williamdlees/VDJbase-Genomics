# write a gene order file, based on the specified bed file

import argparse
import csv
from operator import itemgetter


def main():
    parser = argparse.ArgumentParser(description='Write a gene order file, based on the specified bed files')
    parser.add_argument('bed_file', help='pathname to a bed file listing co-ordinates for the genes')
    parser.add_argument('order_file', help='the gene order file')
    args = parser.parse_args()

    rows = []

    with open(args.bed_file, 'r') as fi:
        reader = csv.DictReader(fi, delimiter='\t', fieldnames=['chain', 'start', 'end', 'gene'])
        for row in reader:
            rows.append({'gene': row['gene'], 'start': int(row['start'])})

    locus_order = sorted(rows, key=itemgetter('start'))
    alpha_order = sorted(rows, key=itemgetter('gene'))

    with open(args.order_file, 'w') as fo:
        fo.write('LOCUS_ORDER = [\n')

        for row in locus_order:
            fo.write('"%s",\n' % row['gene'])

        fo.write(']\n\n')

        fo.write('ALPHA_ORDER = [\n')

        for row in alpha_order:
            fo.write('"%s",\n' % row['gene'])

        fo.write(']\n\n')

        fo.write('PSEUDO_GENES = []')

if __name__ == '__main__':
    main()

