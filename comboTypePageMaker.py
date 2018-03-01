#!/usr/bin/env python
import argparse
import subprocess
import re
import sys
from scipy import special
import numpy

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
                             #choose most common transcript to use a trans. FOR NEW version
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



def hgvsMaker(combo):#makes the HGVSses for all three types

    hgvsses = []
    try:
        cap = re.search('([NMXR]{2}\_\d+\.?\d?)\:c.([-+*]{0,1}\d+.*?)([*\-\>del]{1,4})(.*)',combo)


        hgvsses.append(cap.group(0))
        hgvsses.append(cap.group(1) + ":c." + cap.group(2) + "=" )#WT+
        hgvsses.append(cap.group(1) + ":c.[" + cap.group(2) +  cap.group(3) + cap.group(4) + "];[" + cap.group(2) + "=]")#HET


    except NoneTypeError:

        print "ERROR:HGVS parsing is blown, correct please"

    return hgvsses


def baseHGVS(hgvs):#this is to prevent combining with self. use wT+ as index
    myhgvs = None
    try:
        cap = re.search('([NMXR]{2}\_\d+\.?\d?)\:c.([-+*]{0,1}\d+.*?)([*\-\>del]{1,4})(.*)',combo)
        myhgvs = cap.group(1) + ":c." + cap.group(2) + "=" #WT+

    except NoneTypeError:
        print "ERROR:HGVS parsing is blown, correct please"

    return myhvgs

def dropme(matrix,variant,pos):#this returns true if the var is new, and needs to be juggled.
    isitnew = True

    if variant == matrix[pos]:#
        isitnew = False

    print variant +
    return isitnew

def newMatrix(matx,var,pos):#TRIMS the LINE of the var, keep everything before and delete all after.
    newmat = []                 #trimsthethingy not working correctly

    for oldpos in range(len(matx)):
        if oldpos < pos:
            #newmat[oldpos] = matx[oldpos]
            newmat.append(matx[oldpos])
    return newmat


def comboExpander(variantarray):##maker of the combos, so for all in, produces the uuid and sends to printer
                                #pcant delete the WHOLE variant array. thing, only delete if the variant changes.
    linevals = []               #
    mat =variantarray
    vmat = numpy.mat(variantarray)

    lengthofcombo = len(variantarray)
    lengthofrow = len(variantarray[0])
    numberoftimes = 0

    print vmat.shape
    print special.comb(lengthofcombo * lengthofrow,lengthofcombo,exact = True)

    for row1 in mat[0]:#make this so it deletes only the position that is variable, so, for each, check if it is already in. if not delete the index position that is there
        if len(linevals) < 1:#try init condition
            linevals.append(row1)

        if dropme(mat[0],row1,0):
            linevals = newMatrix(linevals,row1,0)
            linevals.append(row1)


        for row2 in mat[1]:
            if dropme(mat[1],row2,1):
                linevals = newMatrix(linevals,row2,1)
                linevals.append(row2)

            try:#check for third row.

                for row3 in mat[2]:
                    print str(linevals) + "\tEND-LIVE-VALS"
                    print row3 + "\tK WEISEN ROW3 VAL COMING IN NEW"
                    if dropme(mat[2],row3,2):
                        linevals = newMatrix(linevals,row3,2)
                        linevals.append(row3)

                        print str(linevals)  + "\tmat2 and send to printer"

            except IndexError:

                print str(linevals) +  "\tmat1 and send to printer"


#------------main--------##

for combo_ln in muts:
    combo_hgvs_lst = []#contain the types. each row is one variant to use,
    combo_array = combo_ln.split()
    try:
        combos = combo_array[0]
        name = combo_array[1]
        outfi = open(name + '.scf','w')
    except (NameError, IndexError) as e:
        print "ERROR: Input file is wrong"

    #make combos. so, enumerate the types. then, expand them do a triple for each.
    for combination in combos.split(";"):
        combo_hgvs_lst.append(hgvsMaker(combination))

    print combo_hgvs_lst
    comboExpander(combo_hgvs_lst)

























