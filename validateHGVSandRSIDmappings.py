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
vgr = csv.DictReader(open(vgrfi,'r'),delimiter='\t')
muts = open(mutfi,'r')
##-----------defs--------##

def checkMatch(inc_hgvs,inc_chr,inc_pos,matchd):#newtranslated_pos,dictionary

    found = False

    if inc_hgvs in matchd:
        chk_val = "chr" + inc_chr + ":" + inc_pos

        if chk_val == matchd[inc_hgvs]:
            found =True

    return found




#------------main--------##

#blah, turn the hgvs into a dict

lookup = {}

for hg in vgr:
    lookup[hg['hgvs']] = hg['chr'] + ":" + hg['pos']

print("original table loaded.....")


for resln  in muts:

    isitgood = None

    storelocale = {}
    #ok, have to check for valid values, I think the the check for NC might be enough.
    #then, if there is a valid NC value, we can proceed. if nothing at end, then kick out
    #hgvs with a statement of invalid in condition.

    resvalues = resln.split('\t')


    for ncbi in resvalues[1:]:  #keep this, it has NC line parsing.
                                        #sorting is done by version magic.
        if re.match('NC_',ncbi):
            m = re.match('NC_0{1,6}(\d+)\.(\d+):g.(\d+)(.*)',ncbi)
            storelocale[int(m.group(2))] = m.group(3) #loads the version, so that it can be sorted on.
            #print(m.group(3) + "\t" + m.group(2) + "\t" + m.group(1))
            #print("match " + m.group(2) )

            if len(storelocale.keys()) == 1:
                storek =  storelocale.keys()#this statement, will autosort??

                isitgood = checkMatch(resvalues[0],m.group(1),storelocale[storek[0]],lookup)

    if isitgood == True:
        print("MATCH_FOUND\t" + resvalues[0] + "\tchr" + m.group(1) + "\t" + storelocale[storek[0]])

    else:
        print("NOT_FOUND or NOT_RETURNED\t" + resln )





