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
parser.add_argument("--paid", help="the paids to map to",required=True)
parser.add_argument("--gene",help="the NM_ id to use in the file",required=True)
args = parser.parse_args()

mutfi = args.paid
paids = open(mutfi,'r')
toloadfile = args.pgx
loadfile = open(toloadfile,'r')
transcript = args.gene

##-----------defs--------##

def varPrint(row,value):#
    print row + "\t" + value

##--------main----------##


pas=[]

for pa in paids:
    pas.append(pa)

for defs in loadfile:
    cols = defs.split(',')

    hap = "".join(cols[2].split())
    for drug in pas:
        print cols[3] +" ; " + cols[6] +","+ str(cols[7].strip()) +"," + str(drug.strip()) + ","  + hap + "||" + transcript + "," + cols[4] + "," + cols[5]
