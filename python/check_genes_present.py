# Check that all genes in the order file, except pseudogenes, are present in the co-ordinate files

import csv
from importlib.machinery import SourceFileLoader
import os


# files containing the genes in a locus (and pseudogenes to be ignored)
# nb these are files containing python code
order_files = {
    'IGH': 'gene_order/igh_gene_order.py',
}

# the directory containing co-ordinate files
coord_dir = '2021-12-13_franken_gene_coords'

# the loci to check
loci = ['IGH']

# loci containing D genes
heavy_loci = ['IGH', 'TRB']

# co-ordinate files to check for genes
coord_files = {
    'V': ['exon_1.bed', 'exon_2.bed', 'gencode_intron.bed', 'genes.bed', 'heptamer.bed', 'nonamer.bed', 'spacer.bed'],
    'D': ['exon_1.bed', 'genes.bed', 'heptamer.bed', 'nonamer.bed', 'spacer.bed'],
    'J': ['exon_1.bed', 'genes.bed', 'heptamer.bed', 'nonamer.bed', 'spacer.bed'],
}

# enumerate genes in a file
def genes_in_coord_file(filename):
    genes = []
    with open(os.path.join(coord_dir, filename), 'r') as fi:
        reader = csv.reader(fi, delimiter='\t')
        for row in reader:
            genes.append(row[3])

    return genes


# enumerate genes of a type in a locus, excluding pseudogenes
def genes_in_locus(locus):
    gene_order = SourceFileLoader('gene_order', order_files[locus]).load_module()
    genes = list(set(gene_order.LOCUS_ORDER) - set(gene_order.PSEUDO_GENES))

    gene_set = {}
    for gene_type in ['V', 'D', 'J']:
        gene_set[gene_type] = [gene for gene in genes if (locus in gene and gene_type in gene)]

    return gene_set


# run the checks
for locus in loci:
    gene_set = genes_in_locus(locus)

    gene_types = ['V', 'J']

    if locus in heavy_loci:
        gene_types.append('D')

    genes_with_coords = {}
    for gene_type in gene_types:
        for coord_file in coord_files[gene_type]:
            if coord_file not in genes_with_coords:
                genes_with_coords[coord_file] = genes_in_coord_file(coord_file)

    for gene_type in gene_types:
        for gene in gene_set[gene_type]:
            missing = []
            for coord_file in coord_files[gene_type]:
                if gene not in genes_with_coords[coord_file]:
                    missing.append(coord_file)

            if missing:
                print(f"{gene} missing from {', '.join(missing)}")



