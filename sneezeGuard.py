#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import csv
from Bio import SeqIO


parser = argparse.ArgumentParser(description=" match exact Nmer to database of same overlapping nMERs")
parser.add_argument("--qs",help="the file of query, the items you want to check",required=True)
parser.add_argument("--db",help="file of fasta formatted database sequences. Typically known allergens/toxins",required=True)
parser.add_argument("--mer",help="Size of chunk to use for lookup, default 8mer",default=8)


args = parser.parse_args()
qr_fi = args.qs
db_fi = args.db



######++++++++++++++++defs+++++++++++++=#####

def merMaker(size,sequence,seqdct):
    
    offset = 0 

    while offset  <= len(sequence.seq) - args.mer:
        
        seq = str(sequence.seq)

        subseq = seq[int(offset):int(int(offset) + int(size))]
        seqdct[subseq] = sequence.id + " " + sequence.description
           
        offset = offset +1
        




def makeSet(size,fasta):
    seq_dct = {}

    with open(fasta,'r') as fasta_handle:
        for sequence in SeqIO.parse(fasta_handle, "fasta"):
        
           merMaker(size,sequence,seq_dct) 
            

    return(seq_dct)




#####++++++++++++++++++dictionaries++++#####

database = {}#
queries = {}#

###-----------------MAinely-----------#####

database = makeSet(args.mer,db_fi)
queries = makeSet(args.mer,qr_fi)


for ck_mer in queries:
    if ck_mer in database.keys():
        print(ck_mer + "\tFROM\t" + queries[ck_mer]  +"\tIN DATABASE\t" + database[ck_mer])


