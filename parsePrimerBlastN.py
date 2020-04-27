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

def getLongest(qry, prmp,st,sp):#$qryname, primer dict, start and $ stop on query, calculate proportion of match

    length = float((int(sp) - int(st))/prmp[qry])

    return float(length)



def addNew(best,blast,new,pmap):#besthit_dict, blastline if new is True, it will add to the hash, if False, it will just replace the values there.


    if new == True:
        best[blast['queryacc.ver']] = {}
        best[blast['queryacc.ver']][blast['subjectacc.ver']] = {}

    best[blast['queryacc.ver']][blast['subjectacc.ver']] = {"longesthit":getLongest(blast['queryacc.ver'],pmap,blast['q.start'],blast['q.end']), "gaps":int(blast['mismatches'] + blast['gapopens'])}


def queryDict(queryd):

    for key in queryd:
        print(queryd[key])


def Printer(line):

    print(';'.join(line))




def Collapser(bt):

    for res in bt:

        prnstore = []
        prnstore.append(res)

        for qrres in bt[res]:
            prnstore.append( qrres +  ";" +  bt[res][qrres])

        Printer(prnstore)


def findPrimerMatches(bst):#will flip the value, and report the primer hits and what the scores are for subject

    sbjhit = {}

    for qr in bst:
        for subject in bst[qr]:
            sbjhit[subject] = {}
            sbjhit[subject][qr] = str(bst[qr][subject]['longesthit']) + ":" + str(bst[qr][subject]['gaps'])



    Collapser(sbjhit)


##===================MAIN===============##
#get the node into a dict.
pmap = {}#contains the primer lengths.
best = {}#primer is key, value is hit. for every key val, store length of primer hit, and mismatches + gaps.


for primer in pmapfi:
    pkey = primer.split("\t")
    pmap[pkey[0]] = int(pkey[1])



for blast in blncsv:
    #print(blast)

    if blast['queryacc.ver'] in best:

        ##maybe translate subject accesion to species, and dont care about anything less,
        #queryDict(best)


        if blast['subjectacc.ver'] in best[blast['queryacc.ver']]: #query subj value exist

            #print(best[blast['queryacc.ver']] )
            #the alg:
            best_long = best[blast['queryacc.ver']][blast['subjectacc.ver']]['longesthit']
            new_long = getLongest(blast['queryacc.ver'],pmap,blast['q.start'],blast['q.end'])
            best_gap = best[blast['queryacc.ver']][blast['subjectacc.ver']]['gaps']
            new_gap = int(blast['mismatches'] + blast['gapopens'])

            if new_long > best_long and best_gap <= new_gap:
                addNew(best,blast,False,pmap)

                #I may add a condition.




    else:#adding new query, if no query add all in.

        addNew(best,blast,True,pmap)


findPrimerMatches(best)

#print out blast, but before I think I can collect, for each subject, the queries it gets.













