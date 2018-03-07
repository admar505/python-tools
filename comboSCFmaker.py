#!/usr/bin/env python
import argparse
import subprocess
import re
import csv
import numpy
#exporter:https://vapor.veritasgenetics.com/?q=admin/structure/views/view/yt_default_search_solr/edit
#but really:https://vapor.veritasgenetics.com/?q=yt-get-all-variant-links-export&eid=5&return-url=yt-get-all-variant-links-export
parser = argparse.ArgumentParser(description='NOTE:after, convert \'\;\' to comma\n takes the vapor dump to make hte update look up file,\nheader used:Title, transcript, mutation, chrome38, pos38, chrome37, pos37, protein_mutation, uuid')
parser.add_argument("--vap", help="Batch output",required=True)
args = parser.parse_args()

mutfi = args.vap
muts = open(mutfi,'r')
vareader = csv.DictReader(muts,quotechar="\"")
##-----------defs--------##

def uuidPrep(uuid):

    cuttr = re.search('\w+\_\w+\|\|(.*)',uuid)
    return cuttr.group(1)

def cutnid(nid):

    onlyid = re.search('.*?\/(\d+)$',nid)
    return onlyid.group(1)


def varPrint(rows):#HERE:expand to all vartypes, and print out.
                             #(het, homovar, wt+) strategy, it will be key in dict is print tab
    printer = []

    lookups = uuidPrep(rows['uuid'])
    #nid = cutnid(rows['view link'])
    nid = rows['Nid']
    printer.append(lookups)
    printer.append("https://vapor.veritasgenetics.com/?q=node/" + nid)
    printer.append(rows['transcript'] + ":" + rows['Title'])

    print "\t".join(printer)



#------------main--------##

for line in vareader:
    varPrint(line)



































