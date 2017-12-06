#!/usr/bin/env python

import subprocess,sys,os,re,fileinput,argparse,math
import vcf
sys.path.append('/home/nucleo/python-tools/AltObject')
import altobject

parser = argparse.ArgumentParser(description="adds the position info into either the ID col or in INFO")
parser.add_argument("--vcf",help="the vcf file, has to be bgzipped and indexed",required=True)
parser.add_argument("--t",help="the threshold for identifying issues with merged files. if not present or set to zero (default, it will generate training stats)",default=0)


args = parser.parse_args()

vcffi = args.vcf
limit = args.t
#special load for vcffi? (vcffi,"r")

vcf_full = vcf.Reader(open(vcffi,'r'))
#repair for one sample and not two.
tempheader = open('vcf.tmp','w')
headerrewritten = open('tmp.vcf','w')
#subprocess.call(['zgrep \"#\" '  + vcffi],shell=True,stdout=tempheader)
#subprocess.call(['sed -e \'s/unknown\tSample1/Sample1/\' vcf.tmp'],shell=True,stdout=headerrewritten)
#i


#--subs--#
def returnAlts(samples):
    #return None
    my_alt = altobject.AltObj()
    print my_alt.function()
    print my_alt.ab
#--main--#

for var in vcf_full:
    #print var.ALT[0]

    altDict = {}#hold alternates for

    check_for_none = var.ALT[0]

    if check_for_none is not None:
        print var.samples[0]['GT']
            #from here, toss the variant into the alt class, and create alt objects.

        altDict = returnAlts(var)


