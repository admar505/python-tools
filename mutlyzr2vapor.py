#!/usr/bin/env python
import argparse
import subprocess
import re
import csv
import numpy

parser = argparse.ArgumentParser(description='Translate mutalyzer\'s batch output to VAPor input')
parser.add_argument("--mut", help="Batch output",required=True)
parser.add_argument("--tr", help="Good transcripts")
parser.add_argument("--vg", help="the number to start with for VG transcripts, as in 301 or something",type = int)
args = parser.parse_args()

newtrans = open('NEW.vapor.transcripts.txt','w+')
vgTransID = args.vg
mutfi = args.mut
muts = open(mutfi,'r')
transfi = args.tr
trans = open(transfi,'r')
##-----------defs--------##


def varPrint(row1,pos,local):#HERE:expand to all vartypes, and print out.
                             #(het, homovar, wt+)
    np = re.match('(\w\..*?\d+)(\w+)\>(\w+)',pos[1])
    npwt = np.group(1) + np.group(2) + "="
    print row1[0] +","+ pos[1] + "," + pos[0] + "," +   local + "," + pos[1] + "||"+ pos[0]
    print row1[0] +","+ npwt + "," + pos[0] + "," +   local + "," + npwt + "||"+ pos[0]
    print row1[0] +","+ pos[1] + ":0/1," + pos[0] + "," +   local + "," + pos[1] + ":0/1||"+ pos[0]





def bestXX(trarr):#pick closest TR
    toreturn = {}
    for tr in trarr:
        pos = re.match('[cn].[*-]{0,1}(\d+)',tr[1])
        toreturn[int(pos.group(1))] = tr
    keyssorted = sorted(toreturn.keys())
    return toreturn[keyssorted[0]]

def selectAlt(choices,rsid_current,transfile):#array NR > XR > XM

    goodalts = []
    NR = []
    XR = []
    XM = []
    #print str(choices) + "   in select"
    for ind in choices:
        #print vgTransID NEXT:: get the tacos into this space here.
        #ok. lets do this, scan through, take best, and send it back, in best arrray version order.
        if re.match('NR',ind[0]):
            NR.append(ind)
            #print ind
        elif re.match('XR',ind[0]):
            XR.append(ind)
            #print ind
        elif re.match('XM',ind[0]):
            XM.append(ind)
            #print ind

    #now, sort through the stuff.
    if len(NR) > 0:
        goodalts.append(bestXX(NR))
    elif len(XM) > 0:
        goodalts.append(bestXX(XM))
    elif len(XR) > 0:
        goodalts.append(bestXX(XR))
    else:
        goodalts.append(None)
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

def newTrans(vgid,rs,newtrans):#title,trans,version,gene: return array of hgvs for new var.
    newid = "VG_000" + str(vgid)
    newtrans.write('\"' + rs + " " + newid + '\",\"' + newid + '\",\"1\",\"' + rs + '\"\n')
    return newid

def newHGVS(place,vgTransID,rsid, newtrans):
    newvar = []
    newvar.append(newTrans(vgTransID,rsid,newtrans))
    newvar.append('c.10' + place)
    return newvar


def perf(ind):
    return "NULL"

#------------main--------##
goodtrans = {}
for tr in trans:
    m = re.match('(NM_\d+)\.(\d+)',tr)
    goodtrans[m.group(1)] = m.group(2)

for rsid in muts:
    onerow = rsid.split()
    ncbi_ids = onerow[1]
    storelocale = {}#hash, key = version, val = position
    localeprint = ""
    transhgvs = []#store NM ids
    althgvs = [] #store alternate, as in NR,XR ids
    failhgvs = ""#only use if no hgvs is available, ie, only if
    for ncbi in ncbi_ids.split('|'):
        #print ncbi
        if re.match('NC_',ncbi):
            m = re.match('NC_0{1,6}(\d+)\.(\d+):g.(\d+)(.*)',ncbi)
            storelocale[int(m.group(2))] = m.group(3)
            failhgvs = m.group(4)
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
            #perf(1)
            varPrint(onerow,positions,localeprint)

    elif althgvs:
        altID = []
        altID = selectAlt(althgvs,onerow[0],newtrans)
        #print altID
        if altID is not None:
            c = 0
            while c < len(altID):
                varPrint(onerow,altID[c],localeprint)
                c += 1

    else:#NOW: this is working correctly here.)
        varPrint(onerow,newHGVS(failhgvs,vgTransID,onerow[0],newtrans),localeprint)
        vgTransID += 1





