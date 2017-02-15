#!/usr/bin/python

import argparse
import string
import re

parser = argparse.ArgumentParser(description='generate the types for the combos')
parser.add_argument("--haps", help="file with haplotypes to generate.")

args = parser.parse_args()

haps_fi = open(args.haps,"r")


for line in haps_fi.readlines():
    cols = line.split("\t")
    types_in_cols = cols[1].split(",")
    for types in types_in_cols:
        types1 = re.sub("CYP2D6","NM_000106",types)
        types2 = re.sub("CYP2C19","NM_000769",types1)
        types3 = re.sub(";",",",types2)
        print types3 + "\t" + cols[0]
