#!/usr/bin/python
import sys,os,re,fileinput,argparse
import vcf
sys.path.append('/vbin/PyVGRes')
import geneimpacts
import vgr
import csv
from subprocess import call
import pybedtools
import random
#Test if geneImpacts does job I want, and can add to VGR record.
parser = argparse.ArgumentParser(description="VGR and geneimpacts library test!!")
parser.add_argument("--efd",help="fullvcf file",required=True)

args = parser.parse_args()
efdfi = args.efd

vcf_full = vcf.Reader(open(efdfi,'r'))

#-----------------here by DEFSgONS!!-------------------------*





failreturn = vgr.Reader("10")#instanciate empty record

##<---------------++++++++++ MAIN ++++++++++++------------->##
for variant in vcf_full:
    print variant.INFO['ANN']
    failreturn[]
    #newresults = vgr.Reader(open(rsid + '.COMPLETE.txt',"r"))
