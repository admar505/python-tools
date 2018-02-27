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

<<<<<<< HEAD
    print row['uuid']
    parsedvar = re.match('([ATGC]).([ATGC])',row['Title'])
    ref = parsedvar.group(1)
    alt = parsedvar.group(2)
    print ref
    print alt
    print row['Title']
=======
    try:
        parsedvar = re.match('.*?([ATGC]).([ATGC]).*?',row['Title'])
        ref = parsedvar.group(1)
        alt = parsedvar.group(2)

    except AttributeError:
        #parsedvar = re.match('c.(.*)([insdel]+)(.*)',row['Title'])
        parsedvar = re.search('c.(.*?)\:?.*',row['Title'])
        ref = parsedvar.group(1)
        alt = 'delins'

    return row['rsid'] + "-" + ref + "-" + alt
>>>>>>> 4c2743afa45fb5d8c94c54f00f4ac51263348de5

def varPrint(rows):#HERE:expand to all vartypes, and print out.
                             #(het, homovar, wt+) strategy, it will be key in dict is print tab
    printer = {}
                             #THIS is stupid, it comes in threeeeessss
    for var in rows:
        rsid = formatRSID(var)
        #for hgvs in var:
        #print var
        if re.match('.*?:0\/1',var['Title']) is not None:
            printer[9] = "https://vapor.veritasgenetics.com/?q=node/" + var['Nid'].strip()
            printer[8] = var['transcript'] + ":" + var['Title']
        elif re.match('.*=',var['Title']):
            printer[11] = "https://vapor.veritasgenetics.com/?q=node/" + var['Nid'].strip()
            printer[10] = var['transcript'] + ":" + var['Title']
        else:
            printer[7] = "https://vapor.veritasgenetics.com/?q=node/" + var['Nid'].strip()
            printer[6] = var['transcript'] + ":" + var['Title']

    printer[0] = 'chr' + rows[0]['chrome37']

    printer[1] = int(rows[0]['pos37']) - 1
    printer[2] = int(rows[0]['pos37']) + 1
    printer[3] = rsid
    #printer[]


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




































