#!/usr/bin/env python3

import argparse
import string
import re
import sys

parser = argparse.ArgumentParser(description='Removes duplicates from FASTA files.')
parser.add_argument("--fa", help="full fasta file",required=True)

p = parser.parse_args()

fafi = open(p.fa,'r')

fasta = {}

dfs = None   #initialize deflines.

for fln in fafi:

    if ">" in fln:
        dfs = fln.split()[0]
        fasta[dfs] = ""

    else:
        fasta[dfs] = fasta[dfs] + fln

#file.

for line in fasta:

    print(line + "\n" + fasta[line])









