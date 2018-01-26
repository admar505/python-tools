#!/usr/bin/env python

import subprocess,sys,os,re,fileinput,argparse,math
import vcf
sys.path.append('/home/nucleo/python-tools/AltObject')

import loadaltdats

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
tempheader = open('vcf.tmp','w')
headerrewritten = open('tmp.vcf','w')
repeats = open(repeatsfile,'r')

#--subs--#


#--main--#


#for var in vcf_full:
    #print var.ALT[0]
    #print var.CHROM + "\t" + str(var.POS) + "\t" +  str(var.samples[0]['GL'])
    #altDict = {}#hold alternates for

    #check_for_none = var.ALT[0]

    #if check_for_none is not None:
       # print var.samples[0]['GT']
            #from here, toss the variant into the alt class, and create alt objects.



alts = loadaltdats.detRepeats(vcf_full,repeats)

print alts.countRepeats

    #print alts.countRepeats(repeats)
    #for alternatevars in alts:

    #print str(var.POS) + " " + str(alts.assGT_GL) +  " qual " + str(alts.retQUAL)
        #print str(altDict[alternatevars].amivalid)  + "\t" +  str(altDict[alternatevars].getcall) + "\t" +  str(altDict[alternatevars].getindex) +"\t"+  str(altDict[alternatevars].getcall)

