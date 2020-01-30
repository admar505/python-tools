#!/usr/bin/env python3

import argparse
import sys
import gzip
import subprocess
import warnings
sys.path.append('/home/nucleo/biopython')
from Bio import SwissProt

parser = argparse.ArgumentParser(description='given list of UNIPROT ids, get features')
parser.add_argument("--db", help="database dump, gzipped",required=True)
parser.add_argument("--prot", help="desired protein list, should be uniprot ids",required=True)
parser.add_argument("--info",help="information needed from the database, like description as DE, etc",action='append')
args = parser.parse_args()

try:
    swiss = gzip.open(args.db,'rt')
except (TypeError, NameError) as e:
    print('use -h thanks')

prots = open(args.prot,'rt')
features = args.info


#+++++++defs++++++++++++
def recsearch(swissrec,cc):#$$the swiss db record, $$the id desired
    CCfound = None

    for swissrec in rec_db:##OK, SINGLE record.
        for acc in swissrec.accessions:

            if str(cc.strip()) == str(acc):
                CCfound = swissrec

    return CCfound


def desc_dict(full_rec):#DEF creates a dictionary of the stuff from the matching record
                        #description line, super messy though.
    description = {}

    desc_line = full_rec.description.split(";")
    for desc in desc_line:

        try:
            (key,val) = desc.split("=")
            description[str(key.strip())] = val.strip()

        except ValueError:

            description[str(desc.strip())] = desc.strip()

    return description

def printgood(tag_ids,results):

    line = []

    for key in tag_ids:
        line.append(results[key])

    return ','.join(line)

#+++++++main++++++++++++
rec_db = None

try:
    rec_db = list(SwissProt.parse(swiss))
except NameError:
       parser

warnings.warn("\n\n\ndictionary constructed, retrieving feature information.....\n\n")

for prot in prots:

    tag_lst = {} #the final output print holder
    prot_rec = recsearch(rec_db,prot)#matching record, if any
    desc_features = {} #dict of description features

    if prot_rec is not None:#in each feature caller, construct a feature dictionary for the
                            #description line, and make a key value pair with the info,
                            #EC=1.1.1.1  will be dict[EC] = 1.1.1.1
        desc_features = desc_dict(prot_rec)

        for tag in features:

            if tag in desc_features:
                tag_lst[tag] = desc_features[tag]

            else:
                tag_lst[tag] = "NOT_FOUND"

        print_line = printgood(features,tag_lst)

        print(prot.strip() + "," + print_line)


    else:
        returner = ",".join(features)
        print(str(prot.strip()) + "," + returner)



