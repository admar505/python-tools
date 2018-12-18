#!/usr/bin/env python
import sys,os,re,fileinput,argparse
sys.path.append('/home/nucleo/lib/PyVGRes')
import vgr
import csv
import random
parser = argparse.ArgumentParser(description="From a VGR file, make a VAP upload file. need: vgr with VEPHGVS and EFF_PROT, and the ALT and REF present")
parser.add_argument("--vgr",help="vgr file, bgzipped and tabixed",required=True)

args = parser.parse_args()
vcffi = args.vgr

vcf_full = vgr.Reader(open(vcffi,'r'))

#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*
def uuidm(ln):
    vals = []

    vals = ln['VEP_HGVS'].split(":")

    return vals


#####----------------MAIN--------------####      #####----------------MAIN--------------####

for line in vcf_full:

    chrgr = re.search('chr(\w+)',str(line.CHROM))

    print uuidm(line.INFO)[1] + "||" + uuidm(line.INFO)[0] + "," + uuidm(line.INFO)[1] + "," + uuidm(line.INFO)[0] + "," + uuidm(line.INFO)[1] +  ",,," + chrgr.group(1) + "," + str(line.POS) + "," + line.INFO['EFF_PROT']

























