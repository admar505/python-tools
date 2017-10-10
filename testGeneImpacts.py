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
#Plan: take
#
parser = argparse.ArgumentParser(description="VGR and geneimpacts library test!!")
parser.add_argument("--efd",help="snpEff annotated vcf file",required=True)
parser.add_argument("--vep",help="VEP annotated vcf file",required=True)

args = parser.parse_args()
efdfi = args.efd
vepfi = args.vep
efd_full = vcf.Reader(open(efdfi, 'r'))
vep_full = vcf.Reader(open(vepfi, 'r'))
newresults = vgr.Writer(open('test.COMPLETE.txt',"w"))
#-----------------here by DEFSgONS!!-------------------------*
def getGT(sample):#expanded so can be used to extract any sample GT info needed.
    return sample['GT']
def getDP(sample):
    return sample['DP']

def addNewRecord(varrec):#Returns a vgr  record:
    newRes = vgr.model._Record(varrec.CHROM,varrec.POS,varrec.REF,varrec.ALT,{})
    newRes.INFO['QUAL'] = varrec.QUAL
    if 'FB_Genotype' not in newRes.INFO:#one sample will ALWAYS exist 
        newRes.INFO['FB_GenoType'] = getGT(varrec.samples[0])#attempting to get GT, needed big time.
        newRes.INFO['FBTotalDepth'] = getDP(varrec.samples[0])
    if len(varrec.samples) > 1 and 'VS_Genotype' not in newRes.INFO:
        newRes.INFO['VS_GenoType'] = getGT(varrec.samples[1])#
        newRes.INFO['VSTotalDepth'] = getDP(varrec.samples[1])
    #newRes.add_info BLOCK TO HERE APPROVED< GETS GT OK
    return newRes.INFO

    



    
    
    
    
    
    
##<---------------++++++++++ MAIN ++++++++++++------------->##
resData = {}#key is record: val is the VGR line!
for variant in vep_full:
    #resData[str(variant.CHROM) + str(variant.POS) +  str(variant.REF) + str(variant.ALT)] = addNewRecord(variant)
    resData[str(variant.CHROM) + str(variant.POS) +  str(variant.REF) + str(variant.ALT)] = addNewRecord(variant)
    #print variant.INFO['ANN']
    

for row in resData:
    newresults.write_record(resData[row])


#for rec in resData:
#    print rec
#    print resData[rec].INFO['QUAL']
