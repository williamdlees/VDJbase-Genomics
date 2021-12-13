# VDJbase-Genomics

Google doc meeting notes: https://docs.google.com/document/d/1itzTsBQ-0rhJn7LYnXz8YcQFoIKgyS-JrINFz3xalJc/edit

# Structure of github
1. data/: This should contain data that is not expected to change for example the reference
2. results/: This should contain preliminary data/results that we want to share between groups

Each folder should have a date and short description. 

Credit for structure: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424

# 2021-12-08_genomic_gene_sequences
## Intro
I am adding here gene sequences from 157 samples from three cohorts.

## Methods
Sequences from nonomaer to start sites were extracted from the IGenotyper/genomic capture assemblies. 

## Results
The eval_genes.txt file contains 12 columns:
1. Sample name
2. haplotype
3. gene
4. allele
5. RepSeq support
6. heptamer
7. heptamer length
8. spacer
9. spacer length
10. nonomer
11. nonamer length
12. Sequence from nonomer to start site

Example:
```
     1	W-55
     2	h=1
     3	IGHV6-1
     4	IGHV6-1*01
     5	110
     6	CACTGTG
     7	7
     8	CTGGGCTCACACTGACTTCCCCT
     9	23
    10	GGTTTGTGT
    11	9
    12	GTTTGTGTCTGGGCTCACACTGACTTCCCCTCACTGTGTCTCTTGCACAGTAATACACAGCCGTGTCCTCGGGAGTCACAGAGTTCAGCTGCAGGGAGAACTGGTTCTTGGATGTGTCTGGGTTGATGGTTATTCGACTTTTCACAGATACTGCATAATCATTATACCACTTGGACCTGTAGTATGTCCTTCCCAGCCACTCAAGGCCTCTCGATGGGGACTGCCTGATCCAGTTCCAAGCAGCACTGTTGCTAGAGACACTGTCCCCGGAGATGGCACAGGTGAGTGAGAGGGTCTGCGAGGGCTTCACCAGTCCTGGACCTGACTGCTGCAGCTGTACCTGTGACAGGACACCTGGAGACAAAAGGAAACAGCAAAGTGAAACACCCCTCAGTCTGTGAATGCTGCTGTGAATACGGCATCTCCCTGACACTGACCCCATGGGAGGCCCAGCACGGGCAGGAAGATGAGGAAGGAGACAGACA
```

## Conclusions
This is still missing the 5' UTR sequence. I will add it soonish but this is a start. There are some errors in some of the sequences. But those shouldn't matter too much for now since this will be used to build/code the framework.

