#!/usr/bin/env python3


import argparse
import string
import sys
import subprocess
import csv

parser = argparse.ArgumentParser(description='calculate the best primer hit.')
parser.add_argument("--pr", help="primer map, as id '\\t'  length ",required=True)
parser.add_argument("--sp", help='species accession name map, may have to be adjusted for database.', required = True)

parser.add_argument("--bln", help="blastn, parsing the blast hits, from blastn format 6 or 7. add header plz.",required=True)
#header example:
#queryacc.ver   subjectacc.ver  %identity   alignmentlength mismatches  gapopens    q.start q.end   s.start s.end   evalue  bitscore

#approach:
#calculate length of match, store length, and find highest match for the max length
#someday, check for all three, primers + probe matches.
#for every primer, collect the longest match. dont consider evalue, only gaps, and only align percent.
#dont take hit if that condition is not met.

args = parser.parse_args()

blncsv = csv.DictReader(open(args.bln,'r'),delimiter="\t")
pmapfi = open(args.pr,"r")


##+++++++++++++++++++DEFS+++++++++++++++##

def getLongest(pr,st,sp):

    return int(10)



##===================MAIN===============##
#get the node into a dict.
pmap = {}#contains the primer lengths.
best = {}#primer is key, value is hit. for every key val, store length of primer hit, and mismatches + gaps.


for primer in pmap:
    pkey = primer.split("\t")

    pmap[pkey[0]] = pkey[1]




for blast in blncsv:
    print(blast)
    if blast['queryacc.ver'] in best:

        ##translate subject accesion to species, and dont care about anything less,

        if blast['subjectacc.ver'] in best[blast['queryacc.ver']]:

            print(blast['subjectacc.ver'] )


    else:#adding new query, if no query add all in.

        best[blast['queryacc.ver']] = {}
        best[blast['queryacc.ver']][blast['subjectacc.ver']] = {}
        best[blast['queryacc.ver']][blast['subjectacc.ver']] = {"longesthit": getLongest(pmap,blast['q.start'],blast['q.end'])}
        best[blast['queryacc.ver']][blast['subjectacc.ver']] = {"gaps":blast['mismatches'] + blast['gapopens']}





#print out blast













