#!/usr/bin/env python

import subprocess,sys,os,re,fileinput,argparse,math
import vcf
sys.path.append('/home/nucleo/python-tools/AltObject')

from loadaltdats import detGenoType

parser = argparse.ArgumentParser(description="adds the position info into either the ID col or in INFO")
parser.add_argument("--vcf",help="the vcf file, has to be bgzipped and indexed",required=True)
parser.add_argument("--rpt",help="repeat formatted file")
parser.add_argument("--t",help="the threshold for identifying issues with merged files. if not present or set to zero (default, it will generate training stats)",default=0)


args = parser.parse_args()

vcffi = args.vcf
limit = args.t
repeatsfile =args.rpt

#special load for vcffi? (vcffi,"r")

vcf_full = vcf.Reader(open(vcffi,'r'))
#repair for one sample and not two.
#tempheader = open('vcf.tmp','w')
#headerrewritten = open('tmp.vcf','w')
if repeatsfile:
    repeats = open(repeatsfile,'r')

#--subs--#
def homorhet(gt):#if homo, ret true
    homo =  False
    #print gt
    check = gt.split(',')
    if check[0] == check[1]:
        homo = True

    return homo

#--main--#

none = 0
homo = 0
major_alt = 0
minor_alt = 0
extra_alts = 0


for var in vcf_full:
    #print var.ALT[0]
    #print var.CHROM + "\t" + str(var.POS) + "\t" +  str(var.samples[0]['GL'])
    gtglpair = {}
    gtglpair = detGenoType(var)
    how_many_gts = gtglpair.assGT_GL.keys()
    c = 0
    while c < len(how_many_gts):
        gt = how_many_gts[c]
        if homorhet(gt):
            homo += gtglpair.assGT_GL[gt]

        elif c == 1:
            major_alt += gtglpair.assGT_GL[gt]
        elif c == 2:
            major_alt += gtglpair.assGT_GL[gt]
        else:
            extra_alts += gtglpair.assGT_GL[gt]
        c += 1

print homo
print major_alt
print minor_alt
print extra_alt



#alts = loadaltdats.detRepeats(vcf_full,repeats)

#print alts.countRepeats

    #print alts.countRepeats(repeats)
    #for alternatevars in alts:

    #print str(var.POS) + " " + str(alts.assGT_GL) +  " qual " + str(alts.retQUAL)
        #print str(altDict[alternatevars].amivalid)  + "\t" +  str(altDict[alternatevars].getcall) + "\t" +  str(altDict[alternatevars].getindex) +"\t"+  str(altDict[alternatevars].getcall)

