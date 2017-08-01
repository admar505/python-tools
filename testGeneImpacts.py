#!/usr/bin/python
import sys,os,re,fileinput,argparse
import vcf
sys.path.append('/vbin/PyVGRes')
import geneimpacts
import vgr
import csv
from subprocess import call
import random
#Test if geneImpacts does job I want, and can add to VGR record.
parser = argparse.ArgumentParser(description="VGR and geneimpacts library test!!")
parser.add_argument("--efd",help="fullvcf file",required=True)

args = parser.parse_args()
efdfi = args.efd

vcf_full = vcf.Reader(open(efdfi,'r'))

#-----------------here by DEFSgONS!!-------------------------*

def addNewRecord(varrec):
    newRes = vgr.model._Record(varrec.CHROM,varrec.POS,varrec.REF,varrec.ALT,"EFD")
    #newRes.add_info('QUAL',varrec.QUAL)
    #newRes.add_info
    print newRes.INFO
    return newRes
    
    



##<---------------++++++++++ MAIN ++++++++++++------------->##
resData = {}#key is record: val is the VGR line!
for variant in vcf_full:
    resData[str(variant.CHROM) + str(variant.POS) +  str(variant.REF) + str(variant.ALT)] = addNewRecord(variant)
    #print variant.INFO['ANN']
    #newresults = vgr.Reader(open(rsid + '.COMPLETE.txt',"r"))

#for rec in resData:
#    print rec
#    print resData[rec].INFO['QUAL']
