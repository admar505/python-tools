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
targetrs = open(rsidsfi,'r')
trans = open(transfi,'r')
##-----------defs--------##
def nmCheck(ids,trans):
    arr=[]
    ncid = re.match('(NM_\d+)\.\d+:(.*)',ids)
    if ncid.group(1) in trans.keys():
        arr.append(ncid.group(1))
        arr.append(ncid.group(2))
        return arr


#------------main--------##
seentit = {}
rsids = {}
for rsid in targetrs.next():
    rsids[rsid] = rsid

goodtrans = {}
for tr in trans:
    m = re.match('(NM_\d+)\.(\d+)',tr)
    goodtrans[m.group(1)] = m.group(2)

for rsid in muts:
    onerow = rsid.split()
    seentit[onerow[0]] = onerow
    ncbi_ids = onerow[1]
    storelocale = {}#hash, key = version, val = position
    localeprint = ""
    transhgvs = []
    for ncbi in ncbi_ids.split('|'):
        #print ncbi
        if re.match('NC_',ncbi):
            m = re.match('NC_0{1,6}(\d+)\.(\d+):g.(\d+)',ncbi)
            storelocale[int(m.group(2))] = m.group(3)
            if len(storelocale) == 2:
                storek =  storelocale.keys()

                localeprint = m.group(1) +"," + storelocale[storek[0]] +"," + m.group(1) + "," + storelocale[storek[1]]
        elif re.match('NM_',ncbi):
            if nmCheck(ncbi,goodtrans) is not None:
                transhgvs.append(nmCheck(ncbi,goodtrans))
#what if no NMIDS? figure that out::
    for positions in transhgvs:
        print onerow[0] +","+ positions[1] + "," + positions[0] + "," +   localeprint + "," + positions[0] + "||"+ positions[1]

#loser bracket: ensure the rsids are all there.

