#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import vcf
sys.path.append('/home/nucleo/lib/PyVGRes')
import vgr
import csv
from subprocess import call
import pybedtools
import random
#idea: grab the full vcf in a dict/map; grab the omicia vcf in a dict/map;
#then, then make a map to of res. read through the vcf, and add to reslines. If res in vcf, then
#add score and rsid to EFF_EFFECT. if not, then grab line from big vcf, and process like olden times as if new. call recover.RESULTS.txt
#NEW VERSION: go through CSV, pull out the RESults file, if not in, then recreate from vcf.
parser = argparse.ArgumentParser(description="selecting results lines from a list of chr\tpos\tstp fromatted files")
parser.add_argument("--res",help="COMPLETE.txt or whole genome RESULTS.txt file, bgzipped and tabixed as -p vcf",required=True)
parser.add_argument("--skip",help="list of items to skip",required=True)

args = parser.parse_args()
resfi = args.res
skipfi = args.skip

newres = vgr.Writer(open('NEW.' + resfi + '.' + skipfi + ".txt","w"))
res = vgr.Reader(open(resfi,'r'))
skiplist = csv.DictReader(open(skipfi,"r"))


results = {}#stores the results lines;
omiciain = 0
recovered = 0
oneoffed = 0
skipct = 0
#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*

#####----------------MAIN--------------####      #####----------------MAIN--------------####




#convert the skiplist to bed here, allow for readthrough.

skipper = []# pybedtools.bedtool.BedTool()
for skiprow in skiplist:
    start = int(skiprow['pos']) - 1
    stop = int(skiprow['pos']) + 1
    if skiprow['stp'] is not None:
		stop = int(skiprow['stp']) + 1

    skipper.append(skiprow['chr'] + "\t" + str(start)  + "\t" + str(stop))

#convert the omicia csv to annotated bed here:
for region in skipper:
    print region
    region_in_vgr = res.fetch(region.split('\t')[0],int(region.split('\t')[1]),int(region.split('\t')[2]))

    for vgr_record in region_in_vgr:
        newres.write_record(vgr_record)



