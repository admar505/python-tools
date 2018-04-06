#!/usr/bin/env python

import argparse
import numpy as np
from scipy.stats import genextreme
import string
import sys
import subprocess
import random
import resource
import re
import pyfaidx
import numpy

parser = argparse.ArgumentParser(description='HGVS-based Synthetic Read Generator')
parser.add_argument("--mutations", help="file with mutations to generate synthetic reads for")
parser.add_argument("--huref", help="Human Reference FASTA (matching BED) use the faidx (usually alled fai)")
parser.add_argument("--output", help="prefix for file output")
parser.add_argument("--read",help="read length")
parser.add_argument("--dpx",help="read depth")
parser.add_argument("--frag",help="target fragment length")
args = parser.parse_args()
frags = args.frag
mutations_file = args.mutations

# define the reference fasta file
fasta = args.huref

fasta_dct = pyfaidx.Fasta(fasta)#seq is stored as a dict (hash)
#fasta_dct = SeqIO.to_dict(SeqIO.parse(file_fasta,"fasta"))#seq is stored as a dict (hash)

mutations_filehandle = open(mutations_file)
mutations_lines = mutations_filehandle.readlines()
mutations_filehandle.close()

#+++++++defs++++++++++++
def read_name(ID):#provides read name
    readname = "@SYNTH01:NOTREALCCXX:"
    readname += str(random.choice(range(0,12)))
    readname += ":" + str(random.choice(range(0,300))) + ":" + str(random.choice(range(0,1000))) + ":" + str(random.choice(range(0,3000)))

    return readname


def calcCov(chroml,rl,dp):#calc read cov, and number of reads needed.
    dp = int(dp) + 1
    readcount = int(int(chroml)*int(dp)/int(rl))
    return readcount



def genQual(rlen):#generates the qual line
    vals = ['A','F','K','G']
    qual = ''
    for values in range(1,int(rlen)):
        qual  += random.choice(vals)


    return qual

def getSeqs(sq,rl,st):#returns two strings, frag length away from eachother, and facing.
    R1 = sq[st:(int(rl) + int(st))]
    r2start =  int(random.gauss(int(frags),int(rl)/3))
    R2 = ''

    if r2start < len(sq):
        R2 = sq[r2start - int(rl):r2start]

    return [R1,R2]

def createReads(seq,rlen,dpx,idname):#create reads,choose rand selections
    out1 = open(idname + ".R1.fastq",'w')
    out2 = open(idname + ".R2.fastq",'w')

    startsfile = open(idname + "supers",'w')
    term = calcCov(len(seq),rlen,dpx)
    c = 0
    print len(seq)
    while  c < term:
        #startLeft = int(numpy.random.uniform(1,len(seq)))

        if startLeft < len(seq) - int(frags):
            c = c + 1
            rname = read_name(idname)
            reads = getSeqs(seq,rlen,startLeft)
        #   ncheck
            out1.write(rname + "/1\n" + str(reads[0]) + "\n+\n" + genQual(rlen) + "\n")
            out2.write(rname + "/2\n" + str(reads[1]) + "\n+\n" + genQual(rlen) + "\n")


#+++++++defs++++++++++++

for mutations_line in mutations_lines:
    mutations_data = mutations_line.split("\t")
    mchrom = mutations_data[0].strip()
    mleft = int(mutations_data[1].strip())
    moperator = mutations_data[3].strip()
    malteration = mutations_data[4].strip()
    moperator2 = mutations_data[5].strip()
    malteration2 = mutations_data[6].strip()
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
            new_seq = str(old_seq[0:(int(mleft))])+str(old_seq[int(mleft)+(len(malteration)):len(old_seq)])

    if(moperator == ">"):
        new_seq = str(old_seq[0:(mleft - 1)])+malteration+str(old_seq[int(mleft):len(old_seq)])

    if(moperator == "dup"):
        new_seq = str(old_seq[0:mleft])+malteration+str(old_seq[(mleft):len(old_seq)])

    if(moperator == "ins"):
        new_seq = str(old_seq[0:mleft])+malteration+str(old_seq[int(mleft):len(old_seq)])

    if(args.read):
        createReads(new_seq,args.read,args.dpx,new_id)

    else:
        re = re
        new_seq_write = re.sub("(.{50})", "\\1\n", str(new_seq), 0,re.DOTALL)
        OUT.write(">" + mchrom + "\n" + str(new_seq_write) + "\n")

