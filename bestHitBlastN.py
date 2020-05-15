#!/usr/bin/env python3

import argparse
import string
import sys
import subprocess
import csv

parser = argparse.ArgumentParser(description='best reciprocal hit to all the hits, with annotations.')
parser.add_argument("--ann", help="accesion to id lines, so ",required=True)
#parser.add_argument("--sp", help='species accession name map, may have to be adjusted for database.', required = True)
parser.add_argument("--bln", help="blastn, parsing the blast hits, from blastn format 6 or 7. add header plz.",required=True)

#header example:
#queryacc.ver   subjectacc.ver  %identity   alignmentlength mismatches  gapopens    q.start q.end   s.start s.end   evalue  bitscore

#approach:
##calculate length of match, store length, nd find highest match for the max length
##Ithink do legth of match and score?

args = parser.parse_args()

blncsv = csv.DictReader(open(args.bln,'r'),delimiter="\t")
annfi = open(args.ann,"r")


##+++++++++++++++++++DEFS+++++++++++++++##

def getLongest(qry, prmp,st,sp):#$qryname, primer dict, start and $ stop on query, calculate proportion of match

    length = float((int(sp) - int(st))/prmp[qry])
    #print(str(length) + "\t" + str(sp) +"\t" + str(st))
    return float(length)



def addNew(best,blast,newqry,newsbjct,pmap):#besthit_dict, blastline if new is True, it will add to the hash, if False, it will just replace the values there.
                                            #newqry: to add in only new query.
                                            #newsubj:to add in only new sbjct.
    if newqry == True:
        best[blast['queryacc.ver']] = {}


    if newsbjct == True:
        best[blast['queryacc.ver']][blast['subjectacc.ver']] = {}


    #if both subj and qry is already in, just add the new values.
    lhit = getLongest(blast['queryacc.ver'],pmap,blast['q.start'],blast['q.end'])
    gap = int(blast['mismatches']) + int(blast['gapopens'])

    best[blast['queryacc.ver']][blast['subjectacc.ver']] = {"longesthit":lhit, "gaps":gap}

    #print(blast['subjectacc.ver'] + "\t" +str(lhit)  + "\t" + str(gap))


def queryDict(queryd):

    for key in queryd:
        print(key + "+" +  str(queryd[key]))


def Printer(line):

    ct = len(line) - 1
    print(str(ct) + "\t" + ';'.join(line))


def Collapser(bt):

    for res in bt:

        prnstore = []
        prnstore.append(res)

        for qrres in bt[res]:
            prnstore.append( qrres +  ";" +  bt[res][qrres])


        #print(str(res) + "\t"  + str(prnstore))
        Printer(prnstore)


def findPrimerMatches(bst):#will flip the value, and report the primer hits and what the scores are for subject

    sbjhit = {}

    #hmmm, run through, and collect all subjects? ZZ

    for qr in bst:
        #how to flip??#oh, subject needs to go in once, once only. so check and see if it needs to be added.

        for subject in bst[qr]:
            #print(subject + "\t" + qr)

            if subject not in  sbjhit:
                sbjhit[subject] = {}

            sbjhit[subject][qr] = {}
            sbjhit[subject][qr] = str(bst[qr][subject]['longesthit']) + ":" + str(bst[qr][subject]['gaps'])



    #queryDict(sbjhit)
    Collapser(sbjhit)


##===================MAIN===============##
#get the node into a dict.
pmap = {}#contains the annotation headers to sew back on..
bestqr = {}#is key, value is hit. for every key val, store length of primer hit, and mismatches + gaps.
bestdb = {}#primer is key, value is hit. for every key val, store length of primer hit, and mismatches + gaps.


for annline in annfi:
    annkey = primer.split("\s")
    pmap[pkey[0]] = annline

#parse blast file, sort into different buckets.


for blast in blncsv:
    #print(blast)

    if blast['queryacc.ver'] in best:

        ##maybe translate subject accesion to species, and dont care about anything less,
        #queryDict(best)

        #print(best[blast['queryacc.ver']])


        if blast['subjectacc.ver'] in best[blast['queryacc.ver']]: #query subj value exist

            #print(blast['queryacc.ver'] + "\t" + str(best[blast['queryacc.ver']]) )
            #the alg:
            best_long = best[blast['queryacc.ver']][blast['subjectacc.ver']]['longesthit']
            new_long = getLongest(blast['queryacc.ver'],pmap,blast['q.start'],blast['q.end'])
            best_gap = best[blast['queryacc.ver']][blast['subjectacc.ver']]['gaps']
            new_gap = int(blast['mismatches']) + int(blast['gapopens'])
            best_score =

            #print(blast['subjectacc.ver']  +"\t"+ blast['queryacc.ver'] + "\t" + str(new_long) +"\t "+ str(new_gap))


            if new_long >= best_long and best_gap >  new_gap:
                addNew(best,blast,False,False,pmap)#newqry, newsbjct are bools

                #I may add a condition.

        else:#ok, try new logic for adding only the query,
             # I think I might be resetting if I load in new query if I load both
            addNew(best,blast,False,True,pmap)#newqry, newsbjct are bools


    else:#adding new query, if no query add all in.
        addNew(best,blast,True,True,pmap)#newqry, newsbjct are bools


findRecips(best)

#print out blast, but before I think I can collect, for each subject, the queries it gets.


#queryDict(best)










