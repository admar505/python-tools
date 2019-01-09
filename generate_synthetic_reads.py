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
import multiprocessing

parser = argparse.ArgumentParser(description='HGVS-based Synthetic Read Generator')
parser.add_argument("--mutations", help="file with mutations to generate synthetic reads for")
parser.add_argument("--huref", help="Human Reference FASTA (matching BED) use the faidx (usually called fai)")
parser.add_argument("--output", help="prefix for file output")
parser.add_argument("--read",help="read length")
parser.add_argument("--dpx",help="read depth")
parser.add_argument("--frag",help="target fragment length")
parser.add_argument("--threads",help="number of parallel processes to use, default = 2",type=int,default=2)
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
    for values in range(0,int(rlen)):
        qual  += random.choice(vals)

    return qual

def getRevComp(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def getSeqs(sq,rl,st):#returns two strings, frag length away from eachother, and facing.
    R1 = sq[int(st):(int(rl) + int(st))]
    r2start = int(st) + int(random.gauss(int(frags),int(frags)/3))
    R2 = ''

    if r2start < len(sq):
        R2 = sq[r2start - int(rl):r2start]
        R2 = getRevComp(R2.upper())

    return [R1.upper(),R2]


#def runReads(seq,rlen,dpx,idname,out1,out2,term):
#def runReads(seq,rlen,idname,threads,term):
def runReads(prefix):

    #prefix = int(numpy.random.uniform(1,threads))

    out1 = open(idname + "." + str(prefix) +  ".R1.fastq",'w')
    out2 = open(idname + "." + str(prefix) +  ".R2.fastq",'w')

    c = 0
    while c <  loops:
        numpy.random.seed(seed = None)
        startLeft = int(numpy.random.uniform(1,len(seq)))

        if startLeft < len(seq) - int(frags):
            c = c + 1
            rname = read_name(idname)
            reads = getSeqs(seq,rlen,startLeft)

            allnnn1 = re.findall('N',reads[0])
            allnnn2 = re.findall('N',reads[1])

            if len(allnnn1) < len(reads) and len(allnnn2) < len(reads):
                #   ncheck
                out1.write(rname + "/1\n" + str(reads[0]) + "\n+\n" + genQual(rlen) + "\n")
                out2.write(rname + "/2\n" + str(reads[1]) + "\n+\n" + genQual(rlen) + "\n")

    out1.close()
    out2.close()

def createReads(seqs,readlength,dpx,identifiername):#create reads,choose rand selections
    global idname
    idname = identifiername
    global rlen
    rlen = readlength
    global seq
    seq = seqs
    global term
    term = calcCov(len(seq),rlen,dpx)
    #c = 0
    global loops



    loops = int(int(term)/int(args.threads))


    if __name__ == '__main__':
        workers = multiprocessing.Pool(processes=args.threads)#doest seem to

        resultsgetter = workers.map(runReads,range(1,args.threads))

    #FILE handling, make in to one.


    r1 = "full." + str(idname) + ".R1.fastq"
    r2 = "full." + str(idname) + ".R2.fastq"

    totf1 = open(r1,'w')
    totf2 = open(r2,'w')

    subprocess.call("cat " + idname + ".*" +    ".R1.fastq  ",shell=True,stdout = totf1)
    subprocess.call("cat " + idname + ".*" +    ".R2.fastq  ",shell=True,stdout = totf2)
    subprocess.call("rm "   + idname + ".*" +    ".R2.fastq ",shell=True)
    subprocess.call("rm "   + idname + ".*" +    ".R1.fastq ",shell=True)


#+++++++main++++++++++++

for mutations_line in mutations_lines:
    mutations_data = mutations_line.split("\t")
    mchrom = mutations_data[0].strip()
    mleft = int(mutations_data[1].strip())
    moperator = mutations_data[3].strip()
    malteration = mutations_data[4].strip()
    moperator2 = mutations_data[5].strip()
    malteration2 = mutations_data[6].strip()
    old_seq = fasta_dct[mchrom]#whatshould this be? the file handle?
    snpfind = re.search('>',moperator)
    if snpfind is not None:
        mutations_data[3] = 'snp'

    new_id = args.output + "." +  "-".join(mutations_data[0:7]).strip()  +  ".fa"
    if(moperator == "del"):
	#	if (malteration.isalpha()):

        if(moperator2 == "ins"):
            new_seq = str(old_seq[0:(int(mleft))])+str(old_seq[int(mleft)+(len(malteration)):len(old_seq)])#delfirst
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
        OUT = open(new_id,"w")
        new_seq = ""
        new_seq_write = re.sub("(.{50})", "\\1\n", str(new_seq), 0,re.DOTALL)
        OUT.write(">" + mchrom + "\n" + str(new_seq_write) + "\n")

