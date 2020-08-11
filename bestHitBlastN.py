#!/usr/bin/env python3

import argparse
import string
import sys
import subprocess
import csv
from collections import defaultdict
from operator import itemgetter,attrgetter

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

def getLongest(bln):#$qryname, primer dict, start and $ stop on query, calculate proportion of match

    length = int(bln['alignmentlength']) - int(bln['mismatches']) - int(bln['gapopens'])

    return float(length)

def tree():return defaultdict(tree)

def addNew(best,blast,newqry,newsbjct):#edge, it will add
                                            #update: ADD BOTH DIRECTIONS.
                                            #newqry: to add in new first edge
                                            #newsubj:to add in new terminal edge.
                                            #in each, test the opposite.
    #print(str(blast['queryacc.ver']) + "\taddNew line 47")
    lhit = getLongest(blast)

    if newqry == True and newsbjct == False:#newedge add

        if blast['queryacc.ver'] not in best:
            best[blast['queryacc.ver']] = {}

        if blast['subjectacc.ver']  not in best:
            best[blast['subjectacc.ver']] = {}                 #this will be confusing because of the flipped edges, there
                                                                   #is lots of testing are retesting if thinkgs exist.

        #critical error, I have to see things multiple times in order to report, if its once , will never catch it. so I am leaving out the small stuff here.


        best[blast['subjectacc.ver']][blast['queryacc.ver']] = {"longesthit":lhit, "score":float(blast['bitscore'])}
        best[blast['queryacc.ver']][blast['subjectacc.ver']] = {"longesthit":lhit, "score":float(blast['bitscore'])}

    if newsbjct == True and newqry == False:#adding in terminal edge node.

        if blast['queryacc.ver'] not  in best[blast['subjectacc.ver']]:
            best[blast['subjectacc.ver']][blast['queryacc.ver']] = {"longesthit":lhit, "score":float(blast['bitscore'])}

        if blast['subjectacc.ver'] not in best[blast['queryacc.ver']]:
            best[blast['queryacc.ver']][blast['subjectacc.ver']] = {"longesthit":lhit, "score":float(blast['bitscore'])}

        #best[blast['queryacc.ver']][blast['subjectacc.ver']] = {"longesthit":lhit, "score":blast['bitscore']}
        #best[blast['subjectacc.ver']][blast['queryacc.ver']] = {"longesthit":lhit, "score":blast['bitscore']}

    #if both subj and qry is already in, just add the new values.


    if newqry == False and newsbjct == False:#this is if the new vals are better than what is already there.

        best[blast['queryacc.ver']][blast['subjectacc.ver']] = {"longesthit":lhit, "score":float(blast['bitscore'])}
        best[blast['subjectacc.ver']][blast['queryacc.ver']] = {"longesthit":lhit, "score":float(blast['bitscore'])}

    #print(blast['subjectacc.ver'] + "\t" +str(lhit)  + "\t" + str(gap))


def queryDict(queryd):

    for key in queryd:
        print(key + "+DICTIONARY TEST" +  str(queryd[key]))


def Printer(cntg, hit, blastvals,pmp):


    name = pmp[hit]

    scr    = blastvals[cntg][hit]['score']
    length = blastvals[cntg][hit]['longesthit']

    print("FINAL VALUE\t\t" + str(cntg) + "\t" + str(hit) + "\t" + str(scr) + "\t" +  str(length) +"\t" +name)


def Collapser(bt):

    for res in bt:

        prnstore = []
        prnstore.append(res)

        for qrres in bt[res]:
            prnstore.append( qrres +  ";" +  bt[res][qrres])


        #print(str(res) + "\t"  + str(prnstore))
        Printer(prnstore)

def scoreSorted(grph,start):#return best based on score.
                            #how to? get all scores, and sort, and then find the hit that score is too.
    bestedge = ""
    scores = []
    def __sorter__(gr,st,score):

        for hit in gr[st]:

            if float(gr[st][hit]['score']) == score:
                return(hit)

    for hit in grph[start].keys():

        scores.append(float(grph[start][hit]['score']))

    bestedges = sorted(scores,reverse = True)[0]
    besthit = __sorter__(grph,start,bestedges)#for several hours I just missed the "s" at the end.

    return(besthit)


def findRecips(cntg,bst,titles):#cntg is the NODE, bst is edge graph will flip the value, and report the  hits and what the scores are for subject
    print("finding matches")



    #hmmm, run through, and collect all subjects? ZZ

    def __getbestnode4edge__(graph,nde):#graph is the blast list, node is the query
        bestcontig = False

        if nde in graph and bool(graph[nde]) is True:

            bestcontig = scoreSorted(graph,nde)
           # print("     get best node 4 edge " + str(bestcontig))
        return(bestcontig)

    for node  in cntg:
        #print(str(node) + "\tline 111")
        #how to flip??#oh, subject needs to go in once, once only. so check and see if it needs to be added.
        try:

            if node in bst and bool(bst[node]) is True:#filterout the
                bestedge = scoreSorted(bst,node)#get the best match, by values. but which value?
                #print(str(node) +"\t" + str(bst[node])   +"\tBESTEDGE\t" +str(bestedge) )
                bestnode = __getbestnode4edge__(bst,bestedge)#THE BEST TERMINAL NODE!!!

                #print(str(node) + "\tline 112\t" + str(bestedge) +"\t" + str(bestnode))
                if node == bestnode:#THE MAGIC!!! DOES THE NODE match the termimal edges best NODE????
                    Printer(node,bestedge,bst,titles)


        except IndexError:

            print(str(node) + "\tEDGE NOT FOUND\t"  + str(bst[node]))


    #queryDict(bst)
    #Collapser(bjhit)


##===================MAIN===============##
#get the node into a dict.
pmap = {}#contains the annotation headers to sew back on..
bestqr = {}#contig IS SUBJ is key, NODES value is key, STORE ALL CONTIGS, this is NODE list

bestdb = tree()#db and (visa versa  ) EDGES    is key, value is contig/query. for every key val, store length of hit, -  mismatches + gaps.
           #   EDGES!! use as edges and do the best hit versions

for annline in annfi:
    annkey = annline.split(" ")
    pmap[annkey[0]] = annline.strip()

#parse blast file, sort into different buckets.
#so, in short, load, but only take if metrics are better.

for blast in blncsv:
    #print(blast)

    bestqr[blast['subjectacc.ver']] = blast['subjectacc.ver']



    if blast['subjectacc.ver'] in bestdb and blast['queryacc.ver'] in bestdb:#the contig has been seen as key. if not add.

        ##maybe translate subject accesion to species, and dont care about anything less,
        ##ah, so this is going to be harder. can I assume if the contig/query has contains the pair,
        ##then the dbhash does? no, I cannot. I think. just test?

        if blast['subjectacc.ver'] in bestdb[blast['queryacc.ver']] and  blast['queryacc.ver'] in bestdb[blast['subjectacc.ver']]:#query subj value already exist,
        #if blast['queryacc.ver'] in bestdb[blast['subjectacc.ver']]:#query subj value already exist,
                                                                  #so just test if new values are better.
                                                                  #

            best_long = bestdb[blast['queryacc.ver']][blast['subjectacc.ver']]['longesthit']
            #print(str(best_long) + "     BEST LONGEST")
            new_long = getLongest(blast)
            best_score = bestdb[blast['queryacc.ver']][blast['subjectacc.ver']]['score']
            new_score = float(blast['bitscore'])

            #print(blast['subjectacc.ver']  +"\t"+ blast['queryacc.ver'] + "\t" + str(new_long) +"\t "+ str(new_gap))

            score_threshold = best_score*(1-args.better)

            if new_long > best_long and float(new_score) >  float(score_threshold):#THIS PICKS THE BEST,
                addNew(bestdb,blast,False,False)#newqry, newsbjct are bools

                #I may add a condition.


        else:#ok, try new logic for adding only the query,
             # I think I might be resetting if I load in new query if I load both
            addNew(bestdb,blast,False,True)#newqry, newsbjct are bools


    else:#adding new query, if no query add all in.
        addNew(bestdb,blast,True,False)#newqry, newsbjct are bools and just mean

findRecips(bestqr,bestdb,pmap)

#print out blast, but before I think I can collect, for each subject, the queries it gets.
#queryDict(best)










