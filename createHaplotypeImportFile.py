#!/usr/bin/env python
import argparse
import subprocess
import re
import csv
import numpy
#exporter:https://vapor.veritasgenetics.com/?q=admin/structure/views/view/yt_default_search_solr/edit
#but really:https://vapor.veritasgenetics.com/?q=yt-get-all-variant-links-export&eid=5&return-url=yt-get-all-variant-links-export
parser = argparse.ArgumentParser(description='validation of the PGX load. takes the ')
parser.add_argument("--pgx", help="the list from the loading of hte bergermeister input",required=True)
parser.add_argument("--paid", help="the paids to drug name to use",required=True)
parser.add_argument("--lst", help="gene symbol to transcript bed file, usually in NC_PGX/pgx.trans.lst",required=True)

args = parser.parse_args()

mutfi = args.paid
paids = open(mutfi,'r')
toloadfile = args.pgx
loadfile = open(toloadfile,'r')
lstfi = args.lst
translst = open(lstfi,'r')


##-----------defs--------##

def varPrint(row,value):#
    print row + "\t" + value

##--------main----------##


pas={}      #for turning drugs into PAs.
trans={}    #for turning genes into transcript

for pa in paids:
    pcl = pa.split('\t')
    pas[pcl[1].strip()] = pcl[0]


for trn in translst:
    trcl = trn.split('\t')
    trans[trcl[0]] = trcl[1]

for defs in loadfile:
    cols = defs.split(',')
    hap = "".join(cols[2].split())




    try:
        drug = pas[cols[0]]
    except KeyError:
        drug = "NOT_FOUND:" + cols[0]

    try:
        gene = trans[cols[1]]

    except KeyError:
        #print "notfound,likelyaSNP\t" + cols[3] +" ; " + cols[6] +","+ str(cols[7].strip()) +"," + cols[0] + ","  + hap + "||" +cols[1] + "," + cols[4] + "," + cols[5]
        gene = "NOT_FOUND:" + cols[1]

    print cols[3] +" ; " + cols[6] +","+ str(cols[7].strip()) +"," + str(drug.strip()) + ","  + hap + "||" + gene  + "," + cols[4] + "," + cols[5]
