#!/usr/bin/env python
import argparse
import subprocess
import re

parser = argparse.ArgumentParser(description='for building pages from combos, input trans:hgvs;trans:hgvs;trans:hgvs\tlabelID')
parser.add_argument("--mut", help="Batch output",required=True)
#parser.add_argument("--tr", help="Good transcripts")
#parser.add_argument("--vg", help="the number to start with for VG transcripts, as in 301 or something",type = int)
args = parser.parse_args()

#newtrans = open('NEW.vapor.transcripts.txt','w+')
#vgTransID = args.vg
mutfi = args.mut
muts = open(mutfi,'r')
#transfi = args.tr
#trans = open(transfi,'r')
##-----------defs--------##


def varPrint(row1,pos,local):#HERE:expand to all vartypes, and print out.
                             #(het, homovar, wt+)
    #print pos[1]
    #adjust to account for del. Stupid CFTR.

    findsnp = re.compile('(\w\..*?\d+)(\w+)\>(\w+)')
    finddel = re.compile('(\w\..*?\d+)(\w+)del(\w+)')
    findin = re.compile('(\w\..*?\d+)(\w+)del(\w+)')
    global np
    if findin.match(pos[1]) is not None:

        np = findin.search(pos[1])


    elif finddel.match(pos[1]) is not None:

        np = finddel.search(pos[1])


    elif findsnp.match(pos[1]) is not None:

        np = findsnp.search(pos[1])


    #if hgmdfind.search(pos)
    #np = re.match('(\w\..*?\d+)(\w+)\>(\w+)',pos[1])



    npwt = np.group(1) + np.group(2) + "="
    print row1[0] +","+ pos[1] + "," + pos[0] + "," +   local + "," + pos[1] + "||"+ pos[0]
    print row1[0] +","+ npwt + "," + pos[0] + "," +   local + "," + npwt + "||"+ pos[0]
    print row1[0] +","+ pos[1] + ":0/1," + pos[0] + "," +   local + "," + pos[1] + ":0/1||"+ pos[0]


def hgvsMaker(combo):

    hgvsses = []
    for variant in combo.split(";"):
        try:
            cap = re.search('([NMXR]{2}\_\d+\.?\d?\:c.[-+*]{0,1}\d+.*?)([*_\-\>del]){1,4}(.*)',variant)

            print cap.group(1)
            hgvsses.append(cap.group(0))
            hgvsses.append(cap.group(1))

        except NoneType:

            print "ERROR:HGVS parsing is blown, correct please"

    return



def perf(ind):
    return "NULL"

#------------main--------##

for combo_ln in muts:
    combo_hgvs_lst = []
    combo_array = combo_ln.split()
    try:
        combos = combo_array[0]
        name = combo_array[1]
        outfi = open(name + '.scf','w')
    except (NameError, IndexError) as e:
        print "ERROR: Input file is wrong"

    #make combos. so, enumerate the types. then, expand them do a triple for each.

    combo_hgvs_lst.append(hgvsMaker(combos))



