#!/usr/bin/python
import sys,os,re,fileinput,argparse
import vcf
sys.path.append('/vbin/PyVGRes')
sys.path.append('/vbin/goinfo-scrips/AltObject')
import vgr
import altobject
import loadaltdats
import csv
from subprocess import call
#NEW VERSION: go through CSV, pull out the RESults file, if not in, then recreate from vcf.
#chr start   stop    rsid    homohgvs    homourl hethgvs heturl  wthgvs  wturl
parser = argparse.ArgumentParser(description="Assigns to each Type II variant, a vapor URL. if missing, returns error, but keeps going. Also assigns combos. Uses the new RESULTS.txt library!!")
parser.add_argument("--answer",help="the mapping of position to ")
parser.add_argument("--vcf",help="VCF file, bgzipped, indexed with tabix -p vcf")
parser.add_argument("--combo",help="combination types")

args = parser.parse_args()
answerfi = args.answer #-->  to dict, use CSV
vcffi = args.vcf       #--> standard pyvcf
combo = args.combo    # this is array, as this can be several
#special load for vcffi? (vcffi,"r")
newres = vgr.Writer(open('REQUIRED.RESULTS.txt',"w"))#temp name.
resvcf = vcf.Reader(open(vcffi,'r'))
bedfi  = csv.DictReader(open(answerfi,'r'),delimiter='\t')




results = {}#stores the results lines;
#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*





#####----------------MAIN--------------####      #####----------------MAIN--------------####

for bed in bedfi:
    variant = resvcf.fetch(str(bed['chr']),int(bed['start']),int(bed['stop']))
    print variant
    print bed





