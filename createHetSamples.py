#!/usr/bin/python
import sys,os,re,fileinput,argparse
import vcf
sys.path.append('/home/nucleo/lib/PyVGRes')
import vgr
import csv
#idea: grab the full vcf in a dict/map; grab the omicia vcf in a dict/map;
#then, then make a map to of res. read through the vcf, and add to reslines. If res in vcf, then
#add score and rsid to EFF_EFFECT. if not, then grab line from big vcf, and process like olden times as if new. call recover.RESULTS.txt
#NEW VERSION: go through CSV, pull out the RESults file, if not in, then recreate from vcf.
parser = argparse.ArgumentParser(description="replace files with rsid driven nodes and hgvs designations")
parser.add_argument("--res",help="COMPLETE.txt or whole genome RESULTS.txt file")
parser.add_argument("--bed",help="the rsid hgvs url file, using csv.DictReader to parse in")
parser.add_argument("--url",help="the url to use")
parser.add_argument("--hgvs",help="the hgvs to use")
args = parser.parse_args()
resfi = args.res
trnsfi = open(args.bed,'r')
urls = args.url
hg = args.hgvs

res = vgr.Reader(open(resfi,'r'))
#trns = csv.DictReader(open(trnsfi,'r'),delimiter='\t')
newfile = vgr.Writer(open(hg + '.new.' + resfi,"w"))



def rsidline(trans_file,rsid):
    #rets = None
    trans_fi = csv.DictReader(trans_file,delimiter='\t')
    trans_file.seek(0)

    for fline in trans_fi:
        #print fline['rsid'] + "\t" + str(rsid)
        if str(fline['rsid']) == str(rsid):
            #print rsid + "  REPORTED"
            rets = fline

    return rets




for line in res:
    rsold = rsidline(trnsfi,line.INFO['RSID'])
    line.INFO['EFF_HGVS']  = rsold[hg]
    line.INFO['VAPOR_URL'] = rsold[urls]

    newfile.write_record(line)


