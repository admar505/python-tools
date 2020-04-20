#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import csv
import random
parser = argparse.ArgumentParser(description="manipulate HMM records and files")
parser.add_argument("--hmm",help="HMM file",required=True)
parser.add_argument("--chnk",help="Chunk size, number of records per file",required=True)

#number of records per
args = parser.parse_args()
hmmfi = args.hmm
chnk = args.chnk

#psuedo
#
# read through the full file, and count records:
#   everytime the counter hits: and sen it.
#   new file control:
#       cut the hmm off, and a number.
#       increment the file number
#       create a new file when the file reader is done.
#       4717 records. 100 at a time maybe??
#easy peasy.
#

hmm_orig = open(hmmfi,"r")
rct = 0

print(hmmfi)



def outname(filename,currentcount ):

    hmmnaming = hmmfi.split(".")
    currentcount = currentcount + 1

    outfiname = '.'.join(hmmnaming[0:(len(hmmnaming) - 1)]) + "." +  str(currentcount) + ".hmm"

    return outfiname

ficount = 1


outfi = open(outname(hmmfi,rct),"w")


for line in hmm_orig:
    line = line.strip()
    #print(line)

    #check if file needs to be updated
    if rct >= chnk:
        outfi.close()
        outfi = open(outname(hmmfi,rct),"w")



    if line == "//":

        rct = rct + 1
    outfi.write(line)












