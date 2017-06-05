#!/usr/bin/env python
import argparse
import subprocess
import re
import csv

parser = argparse.ArgumentParser(description='Translate mutalyzer\'s batch output to VAPor input')
parser.add_argument("--mut", help="Batch output",required=True)
parser.add_argument("--tr", help="Good transcripts")
parser.add_argument("--rs", help="The rsids you were looking for")
parser.add_argument("--vg", help="the number to start with for VG transcripts, as in 301 or something")
args = parser.parse_args()

vgTransID = args.vg
mutfi = args.mut
muts = open(mutfi,'r')
transfi = args.tr
rsidsfi = args.rs
targetrs = open(rsidsfi,'r')
trans = open(transfi,'r')
##-----------defs--------##
def selectAlt(choices,vgTransID):#array NR > XR > XM
    goodalts = []
    #print str(choices) + "   in select"
    for vrsion in choices:
        if vrsion[0]
        #print vgTransID NEXT:: get the tacos into this space here.
        #ok. lets do this, scan through, take best, and send it back, in best arrray version order.
    return goodalts

def nmCheck(ids,trans):
    arr=[]
    ncid = re.match('(NM_\d+)\.\d+:(.*)',ids)
    if ncid.group(1) in trans.keys():
        arr.append(ncid.group(1))
        arr.append(ncid.group(2))
        return arr


def nrCheck(ids,trans):
    arr=[]
    ncid = re.match('(\w{2}_\d+)\.\d+:(.*)',ids)
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
    transhgvs = []#store NM ids
    althgvs = [] #store alternate, as in NR,XR ids
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
        elif re.match('NR_',ncbi) or re.match('XR_',ncbi) or re.match('XM_',ncbi):
            if nrCheck(ncbi,goodtrans)is not None:
                althgvs.append(nrCheck(ncbi,goodtrans))
##Logic. if no NM, sort other, if none there, use VG_ and create the VG ids.
    #if transhgvs is not None and len(transhgvs) > 0:
    if transhgvs:
        for positions in transhgvs:##This is valid, as thare T>G and T>C or A type mutations.
            #print onerow[0] +","+ positions[1] + "," + positions[0] + "," +   localeprint + "," + positions[0] + "||"+ positions[1]
            print ""

    elif althgvs:
        print str(althgvs) + "\tSTF"
        altID = []
        altID = selectAlt(althgvs,vgTransID)
        print altID
        if altID is not None:
            for alt in altID:
                print altID + ",CATCH NON-NM here"

    else:
        print "NEW TRANS HRERE"
#loser bracket: ensure the rsids are all there.











