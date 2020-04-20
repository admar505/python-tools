#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import csv
import random
parser = argparse.ArgumentParser(description="manipulate HMM records and files")
parser.add_argument("--hmm",help="HMM file",required=True)
parser.add_argument("--chnk",help="Chunk size, number of records per file",required=True,action='append')

#number of records per
args = parser.parse_args()
hmmfi = args.hmm
chunk = args.chnk

#psuedo
#
# read through the full file, and count records:
#   everytime the counter hits: and sen it.
#   new file control:
#       cut the hmm off, and a number.
#       increment the file number
#       create a new file when the file reader is done.
#       4717 records. 100 at a time maybe??
#
#
#
#
#easy peasy.
#



