#!/usr/bin/env python
import sys,os,re,fileinput,argparse
sys.path.append('/home/nucleo/lib/PyVGRes')
import vgr
import csv
import random
parser = argparse.ArgumentParser(description="Generate a lookup table of positions. This was used for the validation of the DV integration. this data will go through the mutalyzer to ensure the position to HGVS matching is correct.")
parser.add_argument("--vgr",help="vgr file, bgzipped and tabixed as -p vcf",required=True)

args = parser.parse_args()
vcffi = args.vgr

vgr_full = vgr.Reader(open(vcffi,'r'))

#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*

def getTrns(full_trans):#create a lookup table of

    returnvals = {}

    for trns in full_trans.split("|"):

        try:
            base,version = trns.split(".")
            returnvals[base] = trns

        except ValueError:
            pass

    return returnvals


def fixHGVS(hg,lookup):

    fixed = None

    try:
        basetrans,cval = hg.split(":")

        if basetrans in lookup:
            fixed = lookup[basetrans] + ":" + cval

    except ValueError:
        pass

    return fixed




#####----------------MAIN--------------####      #####----------------MAIN--------------####

for v in vgr_full:

    if (v.INFO['EFF_HGVS'] or v.INFO['VEP_HGVS']) and v.INFO['EFF_Feature_ID']:

        translist = getTrns(v.INFO['EFF_Feature_ID'])

        #for it was dark. and cloudy all along the land was tacos.

        for term in  v.INFO['EFF_HGVS'].split("|"):
            repairedval = fixHGVS(term,translist)

            if repairedval is not None:#might need a try:katcher here.
                print (v.CHROM + "\t" + v.POS + "\t" + repairedval)






























