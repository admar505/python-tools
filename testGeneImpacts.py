#!/usr/bin/python
import sys,os,re,fileinput,argparse
import vcf
sys.path.append('/vbin/PyVGRes')
import vgr
import csv
from subprocess import call
import pybedtools
import random
#idea: grab the full vcf in a dict/map; grab the omicia vcf in a dict/map;
#then, then make a map to of res. read through the vcf, and add to reslines. If res in vcf, then
#add score and rsid to EFF_EFFECT. if not, then grab line from big vcf, and process like olden times as if new. call recover.RESULTS.txt
#NEW VERSION: go through CSV, pull out the RESults file, if not in, then recreate from vcf.
parser = argparse.ArgumentParser(description="readthrough the omicia CSV, and if needed, reannotate the vcfs for the correct MED.res.files. now ncludes the RESULTS.txt library!!")
parser.add_argument("--csv",help="the omicia csv file(s).")
parser.add_argument("--vcf",help="fullvcf file")
parser.add_argument("--res",help="COMPLETE.txt or whole genome RESULTS.txt file, bgzipped and tabixed as -p vcf")
parser.add_argument("--sample",help="sample name")
parser.add_argument("--skip",help="list of items to skip")
parser.add_argument("--bed",help="the genome bed file to report on")

args = parser.parse_args()
vcf1fi = args.csv
vcffi = args.vcf
resfi = args.res
sample = args.sample
skipfi = args.skip
bedfi = args.bed

#special load for vcffi? (vcffi,"r")
omicia = csv.DictReader(open(vcf1fi,'r'))
vcf_full = vcf.Reader(open(vcffi,'r'))
#vcf_writer = vcf.Writer(open('redo.Merged.vcf', 'w'), vcf_full)
newres = vgr.Writer(open('NEW.' + sample + '.RESULTS.txt',"w"))
res = vgr.Reader(open(resfi,'r'))
reporter = open(sample + ".REPORT.txt","w")
skiplist = csv.DictReader(open(skipfi,"r"))


results = {}#stores the results lines;
omiciain = 0
recovered = 0
oneoffed = 0
skipct = 0
#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*

def LoserRecover(ovcf,rsid):
    failreturn = vgr.Reader("10")#instanciate empty record
    newresults = vgr.Reader(open(rsid + '.COMPLETE.txt',"r"))
    success = 0
    return_value = ""
    for reslines in newresults:
        formattedlines = AddOmicia(ovcf,reslines)
        if formattedlines:
	    formattedlines
            success += 1
            return_value = formattedlines

    if success == 0:
        #print ovcf.ALT
        #failreturn = ovcf.CHROM + ":"  + str(ovcf.POS) +  "\t" +  ovcf.CHROM + "\t" + str(ovcf.POS) + "\t" + str(ovcf.REF) + "\t" + str(ovcf.ALT[0]) + "\tFBRefAlleleCount=0\tFBReferenceAlleleQ=" + str(ovcf.QUAL) + "\tEFF_HGVS=OMICIAUNMAPPABLE:" + ovcf.ID + "\t" + "QUAL=" + str(ovcf.QUAL) + "\t" + "RSID=" + str(ovcf.ID) +  "\n"
        failreturn.CHROM = ovcf[0]
        failreturn.POS = str(ovcf[1])
        failreturn.REF = 'FABRICated'
        failreturn.ALT = 'VAR'
        failreturn.INFO = {}
        failreturn.INFO['VEP_EFFECT'] = "FABRIC_HARD_TO_MAP(" + ovcf[3].rstrip() + ")|(" + ovcf[4].rstrip() + ")"
        failreturn.INFO['FBRefAlleleCount'] = 0
        failreturn.INFO['VEP_HGVS'] = 'UNKNOWN'
        failreturn.INFO['FBReferenceAlleleQ'] = 0
        failreturn.INFO['QUAL'] = ovcf[4].rstrip()
        failreturn.INFO['RSID'] = ovcf[3].rstrip()
        return_value =  failreturn
    return return_value


def AddOmicia(vcf,results):
    if 'VEP_EFFECT' in results.INFO.keys():
        results.INFO['VEP_EFFECT'] = results.INFO['VEP_EFFECT'].rstrip() + "(" + vcf[3].rstrip() + ")|(" + vcf[4].rstrip() + ")"
    elif 'EFF_EFFECT' in results.INFO.keys():
        results.INFO['EFF_EFFECT'] = results.INFO['EFF_EFFECT'] + "(" + vcf[3].rstrip() + ")|(" + vcf[4].rstrip() + ")"
    else:
        results.INFO['VEP_EFFECT'] = "(" + vcf[3].rstrip() + ")|(" + vcf[4].rstrip() + ")"
    results.INFO['RSID'] = vcf[3].rstrip()
    results.INFO['QUAL'] = vcf[4].rstrip()
    return results

def LoserWrite(record,rsid,name):#prot:
    filename = str(rsid) + ".Merged.vcf"
    vcf_writer = vcf.Writer(open(filename, 'w'), vcf_full)
    vcf_writer.write_record(record)
    #print rsid

def LoserReRun(record,rsid,name):
  #  command = "/vbin/GoPipeRUN/goVarAnnotateAndID.sh" + rsid
    #call(["/vbin/GoPipeRUN/goVarAnnotateAndID.sh",rsid])
    effout = open(rsid + '.efd.vcf',"w")
    call(['java -jar /vbin/snpEff/snpEff.jar eff -i vcf -csvStats -hgvs hg19 ' + rsid + '.Merged.vcf'],shell=True,stdout=effout)
