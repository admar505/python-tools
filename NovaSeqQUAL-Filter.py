#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import vcf
import csv
from subprocess import call
import random
#Test if geneImpacts does job I want, and can add to VGR record.
#Plan: take
#
parser = argparse.ArgumentParser(description="Input the vcf, and choose the call qual filter")
parser.add_argument("--vcf",help="snpEff annotated vcf file",required=True)
parser.add_argument("--q",help="VEP annotated vcf file",default=100)
args = parser.parse_args()
vcffi = args.vcf
qthreshold = args.q
vcfname = re.search ('(\S+?).vcf.gz',vcffi)
vcfoutname =  os.path.split(vcfname.group(1))[1] + "-filtered.vcf"
print vcfoutname 

vcfout = open(vcfoutname,'w')
vcf_file = vcf.Reader(open(vcffi,'r'))
outfi = vcf.Writer(vcfout,template=vcf_file)

for var in vcf_file:
    if var.QUAL > qthreshold or var.FILTER is not None:
        outfi.write_record(var) 
        

