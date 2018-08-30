#!/usr/bin/env python
import sys,os,re,fileinput,argparse

parser = argparse.ArgumentParser(description="Prep for VGR processing, this gets rid of NULL\tNULL files")
parser.add_argument("--vgr",help="RESULTS.txt file, VGR formatted",required=True)

args = parser.parse_args()
file = args.vgr
infi = open(file,'r')

offset = 1001

for line in infi:
    lin = re.search('NULL\tNULL',line.strip())
    cols = line.split('\t')
    if lin is not None:

        print cols[0] + "\t" + str(offset) + "\tA\t" + "\t".join(cols[3:len(cols)]).rstrip()
        offset = offset + 1000

    else:
        print line.rstrip()

