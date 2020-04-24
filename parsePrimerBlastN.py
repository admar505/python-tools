#!/usr/bin/env python3


import argparse
import string
import sys
import subprocess
import csv

parser = argparse.ArgumentParser(description='calculate the best primer hit.')
parser.add_argument("--pr", help="primer map, as id '\\t'  length ",required=True)
parser.add_argument("--bln", help="blastn, parsing the blast hits, from blastn format 6 or 7. add header plz.",required=True)
#header example:
#queryacc.ver   subjectacc.ver  %identity   alignmentlength mismatches  gapopens    q.start q.end   s.start s.end   evalue  bitscore

#approach:
#calculate length of match, store length, and find highest match for the max length
#someday, check for all three, primers + probe matches.


args = parser.parse_args()

blncsv = csv.DictReader(open(args.bln,'r'),sep="\t")
pmapfi = open(args.pr,"r")


##+++++++++++++++++++DEFS+++++++++++++++##

##===================MAIN===============##
#get the node into a dict.
pmap = {}#contains the primer lengths.









