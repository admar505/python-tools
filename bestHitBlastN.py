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
parser.add_argument("--better", help="Percentage for score to be better, if within this range it will be replaced if alignment is longer",default=.20,type=float)

#header example:
#queryacc.ver   subjectacc.ver  %identity   alignmentlength mismatches  gapopens    q.start q.end   s.start s.end   evalue  bitscore

#approach:
##calculate length of match, store length, nd find highest match for the max length
##I think do legth of match and score?
##SINCE contigs are what I am attempting to characterize, that
##is what I will start with.



args = parser.parse_args()

blncsv = csv.DictReader(open(args.bln,'r'),delimiter="\t")
annfi = open(args.ann,"r")


##+++++++++++++++++++DEFS+++++++++++++++##

def getLongest(bln,st,sp):#$qryname, primer dict, start and $ stop on query, calculate proportion of match

    length = float((int(sp) - int(st))/[qry])
    #print(str(length) + "\t" + str(sp) +"\t" + str(st))
    return float(length)



def addNew(best,blast,newqry,newsbjct,qry):#besthit_dict, blastline if new is True, it will add to the hash, if False, it will just replace the values there. this is contig blast agnostic, pass it and it will add.
                                            #newqry: to add in only new query.
                                            #newsubj:to add in only new sbjct.
                                            #need to adjust to know which to add, the query or subj as key.
    print(blast)
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


def findRecips(bst):#will flip the value, and report the  hits and what the scores are for subject

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
bestqr = {}#CONTIG is is key, value is hit. for every key val, store length of  hit,score, and mismatches + gaps.
           #for length, subtract gaps and mismatches, so, make it a qualified hit.
           #may need to adjust to allow for gaps. ways to gap but break? alse perhaps HSP merging?
bestdb = {}#db is key, value is contig/query. for every key val, store length of hit, -  mismatches + gaps.


for annline in annfi:
    annkey = annline.split("\s")
    pmap[annkey[0]] = annline

#parse blast file, sort into different buckets.
#so, in short, load, but only take if metrics are better.


for blast in blncsv:
    #print(blast)

    if blast['queryacc.ver'] in bestqr:#the contig has been seen as key. if not add.

        ##maybe translate subject accesion to species, and dont care about anything less,
        ##ah, so this is going to be harder. can I assume if the contig/query has contains the pair,
        ##then the dbhash does? no, I cannot. I think. just test?

        #print(best[blast['queryacc.ver']])


        if blast['subjectacc.ver'] in bestdb[blast['queryacc.ver']] and b:#query subj value already exist,
                                                                  #so just test if new values are better.
                                                                  #
            #print(blast['queryacc.ver'] + "\t" + str(best[blast['queryacc.ver']]) )
            #the alg:
            best_long = bestdb[blast['queryacc.ver']][blast['subjectacc.ver']]['longesthit']
            new_long = getLongest(blast['queryacc.ver'],pmap,blast['q.start'],blast['q.end'])
            best_gap = best[blast['queryacc.ver']][blast['subjectacc.ver']]['gaps']
            new_gap = int(blast['mismatches']) + int(blast['gapopens'])
            best_score = blast['bitscore']

            #print(blast['subjectacc.ver']  +"\t"+ blast['queryacc.ver'] + "\t" + str(new_long) +"\t "+ str(new_gap))


            if new_long >= best_long and best_gap >  new_gap:
                addNew(best,blast,False,False)#newqry, newsbjct are bools

                #I may add a condition.

        else:#ok, try new logic for adding only the query,
             # I think I might be resetting if I load in new query if I load both
            addNew(best,blast,False,True)#newqry, newsbjct are bools


    else:#adding new query, if no query add all in.
        addNew(bestqr,blast,True,True)#newqry, newsbjct are bools


findRecips(bestqr,bestdb,pmap)

#print out blast, but before I think I can collect, for each subject, the queries it gets.


#queryDict(best)










