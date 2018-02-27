#!/usr/bin/env python
import argparse
import subprocess
import re
import csv
import numpy

parser = argparse.ArgumentParser(description='newest version of input, example:rsid,Title,transcript,chrome37,pos37,chrome38,pos38,uuid,Nid\nrs478304,c.10G>T,rVG_000200,11,65494260,11,65726789,c.10G>T||VG_000200,906184')
parser.add_argument("--mut", help="Batch output",required=True)
args = parser.parse_args()

mutfi = args.mut
muts = open(mutfi,'r')
vareader = csv.DictReader(muts)
##-----------defs--------##

def formatRSID(row):

    print row['uuid']
    parsedvar = re.match('([ATGC]).([ATGC])',row['Title'])
    ref = parsedvar.group(1)
    alt = parsedvar.group(2)
    print ref
    print alt
    print row['Title']

def varPrint(rows):#HERE:expand to all vartypes, and print out.
                             #(het, homovar, wt+) strategy, it will be key in dict is print tab
    printer = {}

    rsid = formatRSID(rows)
    for hgvs in rows:
        if re.match('.*?:0\/1',hgvs['Title']) is not None:
          print hgvs['uuid'].strip()

    printer[0] = 'chr' + rows[0]['chrome37']

    printer[1] = int(rows[0]['pos37']) - 1
    printer[2] = int(rows[0]['pos37']) + 1
    printer2 = []

    for val in printer.keys():

        printer2.append(str(printer[val]).strip())



    print '\t'.join(printer2).strip()


#------------main--------##
rsidpack = []#fill with three rsid lines, and push to printer.
for line in vareader:
    if len(rsidpack)== 3:#reset array.
        varPrint(rsidpack)
        rsidpack = []
        rsidpack.append(line)
    else:
        rsidpack.append(line)




































