#!/usr/bin/env python
import argparse
import subprocess
import re
import csv

parser = argparse.ArgumentParser(description='translate mutalyzer\'s batch output to VAPor input')
parser.add_argument("--mut", help="batch output",required=True)
parser.add_argument("--tr", help="good transcripts")
parser.add_argument("--rs", help="the rsids you were looking for")
args = parser.parse_args()

mutfi = args.mut
muts = open(mutfi,'r')
transfi = args.tr
rsidsfi = args.rs
targetrs = csv.DictReader(open(rsidsfi,'r'))
trans = csv.DictReader(open(transfi,'r'))
##-----------defs--------##
#def :


#------------main--------##







