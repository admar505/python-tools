#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import vcf
sys.path.append('/home/diegoisi/PyVGRes')
import vgr
import csv
from subprocess import call
import pybedtools
import random
#idea: make sure all the NMID are what we expect it to be. 
parser = argparse.ArgumentParser(description="NM_0001321.2 type of line, did get it?")
parser.add_argument("--res",help="COMPLETE.txt or whole genome RESULTS.txt file, bgzipped and tabixed as -p vcf",required=True)

args = parser.parse_args()
resfi = args.res

res = open(resfi,'r')

#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*

#####----------------MAIN--------------####      #####----------------MAIN--------------####


#convert the omicia csv to annotated bed here:
for region in res:
    to_find = None 
    first_capture = None
    second_capture = None
    reg = region.split()
    for val in reg:
        NM_ID = re.search('.*?(NM_\d+).*?',val)

        if NM_ID is not None and to_find is None:
            #print NM_ID.group(1)
            to_find = NM_ID.group(1)

        elif NM_ID is not None and first_capture is None:
            first_capture = NM_ID.group(1)

        elif NM_ID is not None and first_capture is not None:
            second_capture = NM_ID.group(1)
    print str(to_find) + "\t" + str(first_capture) + "\t" + str(second_capture)
    if to_find != first_capture and to_find != second_capture:
        print str(to_find) + "\t" + str(first_capture) + "\t" + str(second_capture) + "\t" + region 
    	



