#!/usr/bin/env python
import argparse
import subprocess
import re
import csv
import numpy

parser = argparse.ArgumentParser(description='given a list of hgvss\'/rsids from mutalyzer, see if they the same as the list from the genome  ')
parser.add_argument("--mut", help="Batch output",required=True)
parser.add_argument("--vg", help="the initial file, chr\tpos\tHGVS",required=True)
parser.add_argument("--rsid", help="if the batch (--mut) is rsids, use this flag.",default = False)
args = parser.parse_args()

newtrans = open('NEW.vapor.transcripts.txt','w+')
mutfi = args.mut
vgrfi = args.vg
muts = csv.Reader(open(vgrfi,'r'),delimiter='\t')
##-----------defs--------##



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


def RSIDadder(hgv):

    m = re.search('(\w+)([>delinsdup]{1,3})(\w+)',hgv)#absolutely RETARDED, captures min not max

    if m.group(2) == 'l':
        rsidadder = "-" + "del-" + m.group(3)

    elif m.group(2) == 's':
        rsidadder = "-" + "ins-" + m.group(3)

    elif m.group(2) == 'u':
        rsidadder = "-" + "dup-" + m.group(3)

    else:
        rsidadder = "-" + m.group(1) + "-" + m.group(3)

    return rsidadder




def perf(ind):
    return "NULL"

#------------main--------##

#blah, turn the hgvs into a dict

lookup = {}

for hg in vgr:
    lookup[hg['hgvs']] = hg['chr'] + ":" + hg['pos']




for resln  in muts:

    #onerow = rsid.split()
    #ncbi_ids = onerow[1]
    storelocale = {}#hash, key = version, val = position
    localeprint = ""
    transhgvs = []#store NM ids
    althgvs = [] #store alternate, as in NR,XR ids
    failhgvs = ""#only use if no hgvs is available, ie, only if

    for ncbi in ncbi_ids.split('|'): #keep this, it has NC line parsing.
        #print ncbi
        if re.match('NC_',ncbi):
            m = re.match('NC_0{1,6}(\d+)\.(\d+):g.(\d+)(.*)',ncbi)
            storelocale[int(m.group(2))] = m.group(3)
            failhgvs = m.group(4)

            if len(storelocale) == 2:
                storek =  storelocale.keys()

                localeprint = m.group(1) +"," + storelocale[storek[0]] +"," + m.group(1) + "," + storelocale[storek[1]]







