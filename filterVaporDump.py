#!/usr/bin/env python
import argparse
import subprocess
import re
import csv
import numpy

parser = argparse.ArgumentParser(description='newest version of input, example:rsid,Title,transcript,chrome37,pos37,chrome38,pos38,uuid,Nid\nrs478304,c.10G>T,rVG_000200,11,65494260,11,65726789,c.10G>T||VG_000200,906184')
parser.add_argument("--mut", help="Batch output of vapor",required=True)
parser.add_argument("--lst", help="List of stuff needed, trnx,c.VAL is header",required=True)
args = parser.parse_args()

mutfi = args.mut
lstfi = args.lst
muts = open(mutfi,'r')
vareader = csv.DictReader(muts,quotechar='\"')
lstreader = csv.DictReader(open(lstfi))
outfi = open('vapor.cut.csv','w')
output = csv.DictWriter(outfi,fieldnames=["Title","transcript","uuid","Nid","rsID","PGX Reference CMO","ClinicalSignificanceRank","Phenotype Title","chrome37","chrome38","pos37","pos38"])
##-----------defs--------##

def formatHGVS(row):

    try:
        parsedvar = re.match('(c.[1234567890*\-+ATGC_]+)[=>delins]{1,4}.*',row)
        ref = parsedvar.group(1)

    except AttributeError:
        #parsedvar = re.match('c.(.*)([insdel]+)(.*)',row['Title'])
        print "cant do this:" + str(row)
        ref = None

    return ref

def varPrint(rows):#HERE:expand to all vartypes, and print out.

    return True


#------------main--------##
finder = {}
for find in lstreader:
    finder[formatHGVS(find['c.VAL'])] = find["trnx"]

output.writeheader()

for line in vareader:
    hgvscheck = formatHGVS(line['Title'])
    if hgvscheck in finder:
        if finder[hgvscheck] == line['transcript']:
           #output.writerow(line[["Title","transcript","uuid","Nid","rsID","chrome37","chrome38","pos37","pos38"]])
           #output.writerow(line[{"Title","transcript","uuid","Nid","rsID","chrome37","chrome38","pos37","pos38"}])
           output.writerow(line)

































