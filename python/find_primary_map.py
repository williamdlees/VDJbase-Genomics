# for each read in a sam file extract, if it's a secondary read, find where its primary is located in the whole file

master_sam_file = '../sc-11_ccs_to_ref_phased.sorted.sam'
extract_sam_file = 'visibleData_RG0.sam'
outfile = 'visibleData_RG0_sec_to_primary.csv'

from receptor_utils import simple_bio_seq as simple
import csv

primaries = {}

with open(master_sam_file, 'r') as fi:
    for row in csv.reader(fi, delimiter='\t'):
        if len(row) > 15 and int(row[1]) < 256:
            primaries[row[0]] = int(row[3])

results = []

with open(extract_sam_file, 'r') as fi:
    for row in csv.reader(fi, delimiter='\t'):
        if len(row) > 15 and int(row[1]) >= 256:
            if row[0] in primaries:
                results.append({'id': row[0], 'loc': primaries[row[0]]})
            else:
                results.append({'id': row[0], 'loc': -1})

simple.write_csv(outfile, results)