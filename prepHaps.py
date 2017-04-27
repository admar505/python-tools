
import argparse
import numpy as np
from scipy.stats import genextreme
import string
import sys
import subprocess
import random
import re

parser = argparse.ArgumentParser(description='MErge the haplotype version 2 files, to give lookup table stuff.')
parser.add_argument("--node", help="this contains the node ==> group rep identifier ")
parser.add_argument("--ids", help="this is the representative ==> child haplotypes")
args = parser.parse_args()

node_fi = args.node
ids_fi = args.ids
##+++++++++++++++++++DEFS+++++++++++++++##

##===================MAIN===============##
#get the node into a dict.
parent2node = {}
nodesFH = open(node_fi,"r")
nodes =  nodesFH.readlines()
for node in nodes:
    fields = node.split(',')
    parent2node[fields[0]] = ([fields[1],fields[2].strip()])
nodesFH.close()

hapsFH = open(ids_fi,"r")
haps = hapsFH.readlines()
for hapt in haps:
    cols = hapt.strip().split('\t')
    indhap = cols[2].split(',')
    for hap in indhap:
        vaporid = re.search('\S+\|\|(\w+)',cols[0])
        print  hap + "||" + vaporid.group(1) + "," + parent2node[cols[0]][0] + "," + parent2node[cols[0]][1]

