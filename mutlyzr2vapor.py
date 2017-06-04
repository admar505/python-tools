#!/usr/bin/env python
import argparse
import subprocess
import re
import csv

parser = argparse.ArgumentParser(description='translate mutalyzer\'s batch output to VAPor input')
parser.add_argument("--mut", help="batch output",required=True)
parser.add_argument("--tr", help="good transcripts")
parser.add_argument("--rs", help="the rsids you were looking for")
args = parser.parse_args()

mutfi = args.mut
muts = open(mutfi,'r')
transfi = args.tr
rsidsfi = args.rs
targetrs = csv.DictReader(open(rsidsfi,'r'))
trans = csv.DictReader(open(transfi,'r'))
##-----------defs--------##
#de:


#------------main--------##
seentit = {}
rsids = {}
for rsid in targetrs.next():
    rsids[rsid] = rsid

goodtrans = {}
for tr in trans.next():
    goodtrans[tr] = tr

for rsid in muts:
    
    print rsid.strip()
    onerow = rsid.split()
    seentit[onerow[0]] = onerow
    ncbi_ids = onerow[1]
    storelocale = {}
    localeprint = ""
    for ncbi in ncbi_ids.split('|'):   
        #print ncbi
        if re.match('NC_',ncbi):
            m = re.match('NC_0{1,6}(\d+)\.(\d+):g.(\d+)',ncbi)
            print m.group(2)
            print ncbi
            storelocale[int(m.group(2))] = m.group(3)
            if len(storelocale) == 2:
                storekeys =  storelocale.keys()
                
                print sorted(storekeys)
