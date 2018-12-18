#!/usr/bin/env python
import argparse
import re
from collections import defaultdict

#exporter:https://vapor.veritasgenetics.com/?q=admin/structure/views/view/yt_default_search_solr/edit
#but really:https://vapor.veritasgenetics.com/?q=yt-get-all-variant-links-export&eid=5&return-url=yt-get-all-variant-links-export
parser = argparse.ArgumentParser(description='expands the CYP2D6 uuids to many drugs by using metabolizer status as primary key ')
parser.add_argument("--paid", help="the list from EXEC that has the metabolizer to recommendation list ",required=True)
parser.add_argument("--pgx", help="the uuid to metabolizer list",required=True)

args = parser.parse_args()

mutfi = args.paid
paids = open(mutfi,'r')
toloadfile = args.pgx
loadfile = open(toloadfile,'r')


##-----------defs--------##

def varPrint(row,value):#
    print row + "\t" + value

def tree(): return defaultdict(tree)


##--------main----------##


pas=tree()      #for turning drugs into PAs.

for pa in paids:
    pcl = pa.split('\t')

    pas[pcl[0].strip()][pcl[1].strip()] = pcl



for defs in loadfile:
    cols = defs.split('\t')
    hapstatus = str(cols[1].strip())

    for drug in pas[hapstatus]:

        rec = pas[hapstatus][drug]


        print rec[0].title() + " ; "  + rec[8] + "," + rec[9].strip() + "," + rec[1].strip() + "," + cols[0] + ","  + rec[6] + "," + rec[7]
