# for each read in a sam file extract, if it's a primary read, find where its secondaries are located in the whole file

master_sam_file = '../sc-11_ccs_to_ref_phased.sorted.sam'
extract_sam_file = 'visibleData_RG2.sam'
outfile = 'visibleData_RG2_primary_to_sec.csv'

from receptor_utils import simple_bio_seq as simple
import csv

secondaries = {}

with open(master_sam_file, 'r') as fi:
    for row in csv.reader(fi, delimiter='\t'):
        if len(row) > 15 and int(row[1]) >= 256:
            if row[0] in secondaries:
                secondaries[row[0]].append(str(row[3]))
            else:
                secondaries[row[0]] = [str(row[3])]

results = []

with open(extract_sam_file, 'r') as fi:
    for row in csv.reader(fi, delimiter='\t'):
        if len(row) > 15 and int(row[1]) < 256:
            if row[0] in secondaries:
                results.append({'id': row[0], 'loc': ', '.join(secondaries[row[0]])})
            else:
                results.append({'id': row[0], 'loc': -1})

simple.write_csv(outfile, results)