#!/usr/bin/env python
import argparse
import subprocess
import re
import csv
import numpy
#exporter:https://vapor.veritasgenetics.com/?q=admin/structure/views/view/yt_default_search_solr/edit
#but really:https://vapor.veritasgenetics.com/?q=yt-get-all-variant-links-export&eid=5&return-url=yt-get-all-variant-links-export
parser = argparse.ArgumentParser(description='NOTE:takes the TA repeat and produces the *1 etc type translation')
parser.add_argument("--vap", help="Output",required=True)
args = parser.parse_args()

mutfi = args.vap
muts = open(mutfi,'r')
##-----------defs--------##

def varPrint(rows):#HERE:
    return None




#------------main--------##
star={}

for line in muts:

    (startype,counttype) = line.split()
    star[counttype] = startype

for c in star:
    for d in star:

        st = []
        rpt = []
        rpt.append(c)
        rpt.append(d)
        st.append(star[c])
        st.append(star[d])
        rpt =  sorted(rpt)
        st = sorted(st)

        print "c.-53_-52TA[" + rpt[0] + "]/[" + rpt[1] + "]||NM_000463\t" + st[0] + "/" + st[1] + "||NM_000463"





