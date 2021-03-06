#!/usr/bin/env python
import argparse
import numpy as np
from scipy.stats import genextreme
import string
import sys
import subprocess
import random
import resource
import unittest
from time import sleep
from Bio import SeqIO
import re
import pyfaidx

parser = argparse.ArgumentParser(description='HGVS-based Synthetic Read Generator')
parser.add_argument("--mutations", help="file with mutations to generate synthetic reads for")
parser.add_argument("--huref", help="Human Reference FASTA (matching BED) use the faidx (usually alled fai)")
parser.add_argument("--output", help="prefix for file output")
args = parser.parse_args()

mutations_file = args.mutations


# define the reference fasta file
fasta = args.huref

fasta_dct = pyfaidx.Fasta(fasta)#seq is stored as a dict (hash)
#fasta_dct = SeqIO.to_dict(SeqIO.parse(file_fasta,"fasta"))#seq is stored as a dict (hash)

mutations_filehandle = open(mutations_file)
mutations_lines = mutations_filehandle.readlines()
mutations_filehandle.close()

for mutations_line in mutations_lines:
	mutations_data = mutations_line.split("\t")
	#mchrom = "chr" + mutations_data[1].strip()
	mchrom = mutations_data[0].strip()
       # print(mchrom)
	mleft = int(mutations_data[1].strip())
	#mright = int(mutations_data[2].strip())
	moperator = mutations_data[3].strip()
	malteration = mutations_data[4].strip()
	moperator2 = mutations_data[5].strip()
	malteration2 = mutations_data[6].strip()
        #old_seq = fasta_dct[mchrom].seq#whatshould this be? the file handle?
        old_seq = fasta_dct[mchrom]#whatshould this be? the file handle?
        new_id = args.output + "." +  "-".join(mutations_data[0:7]).strip()  +  ".fa"
        OUT = open(new_id,"w")
	new_seq = ""
	if(moperator == "del"):
	#	if (malteration.isalpha()):
		if(moperator2 == "ins"):
                    new_seq = str(old_seq[0:(int(mleft))])+str(old_seq[int(mleft)+(len(malteration)):len(old_seq)])
                    new_seq = str(new_seq[0:mleft])+malteration2+str(new_seq[int(mleft):len(old_seq)])
                else:
        #           print("OOHH NOOOEESS")
                    new_seq = str(old_seq[0:(int(mleft))])+str(old_seq[int(mleft)+(len(malteration)):len(old_seq)])
                #print("INDEED . TREEEEUUUOOEEE")
                #print(new_seq)

	if(moperator == ">"):
                new_seq = str(old_seq[0:(mleft - 1)])+malteration+str(old_seq[int(mleft):len(old_seq)])
	if(moperator == "dup"):
#		new_seq = old_seq[0:rel_pos-1]+malteration+old_seq[(rel_pos-1):len(old_seq)]
                new_seq = str(old_seq[0:mleft])+malteration+str(old_seq[(mleft):len(old_seq)])
                #print("DUP")
	if(moperator == "ins"):
#		new_seq = old_seq[0:rel_pos-1]+malteration+old_seq[(rel_pos-1):len(old_seq)]
                new_seq = str(old_seq[0:mleft])+malteration+str(old_seq[int(mleft):len(old_seq)])
	#print(str(len(new_seq)))
	#print(new_seq)
	#print(str(len(old_seq)))
	#print(old_seq)
#	#amp_seq[NAME] = new_seq
        #print(new_seq)
        re = re
        new_seq_write = re.sub("(.{50})", "\\1\n", str(new_seq), 0,re.DOTALL)
        OUT.write(">" + mchrom + "\n" + str(new_seq_write) + "\n")
        #SeqIO.write(str(new_seq),OUT,"fasta")


