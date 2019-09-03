#!/usr/bin/python
import sys,os,re,fileinput,argparse
import vcf
sys.path.append('/home/diegoisi/lib/PyVGRes')
import vgr
import csv
#idea: grab the full vcf in a dict/map; grab the omicia vcf in a dict/map;
#then, then make a map to of res. read through the vcf, and add to reslines. If res in vcf, then
#add score and rsid to EFF_EFFECT. if not, then grab line from big vcf, and process like olden times as if new. call recover.RESULTS.txt
#NEW VERSION: go through CSV, pull out the RESults file, if not in, then recreate from vcf.
parser = argparse.ArgumentParser(description="replace files with rsid driven nodes and hgvs designations, if the variant DOESNT have rsid, then remove from file thanks")
parser.add_argument("--res",help="COMPLETE.txt or whole genome RESULTS.txt file")
parser.add_argument("--bed",help="the answer.bed, but with the rsid-REF-ALT converted to rsid")
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

def correctGT(req_gt,bedline):
    genotype = None

    def __getALT__(rsid,allele):

        alt = None
        fgt = re.search('rs\w+\-(\w+)\-(\w+)',rsid)

        if "hom" in allele:
            alt = fgt.group(2) + fgt.group(2)

        elif "wt" in allele:
            alt = fgt.group(1) + fgt.group(1)

        else:
            alt = ''.join(sorted(''.join([fgt.group(1),fgt.group(2)])))

        return(alt)
    #------#
    genotype = __getALT__(bedline['rsid'],req_gt)

    return genotype



def rsidline(trans_file,rsid):#returns the correct var line, but, I have the special line.
    rets = None
    trans_fi = csv.DictReader(trans_file,delimiter='\t')
    trans_file.seek(0)

    for fline in trans_fi:
        #if str(fline['rsid']) == str(rsid):
        if str(rsid) in str(fline['rsid']):
            rets = fline

    return rets


for line in res:
    if line.ALT != "NULL":
        rsold = rsidline(trnsfi,line.INFO['RSID'])
        line.INFO['EFF_HGVS']  = rsold[hg]
        line.INFO['VAPOR_URL'] = rsold[urls]
        line.INFO['FBGenoType'] = correctGT(hg,rsold)

        newfile.write_record(line)

    else:

        newfile.write_record(line)


