#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import csv
import random
parser = argparse.ArgumentParser(description="takes a file of CYC data, and produces pairwise info for cytoscape network viewing" )
parser.add_argument("--fi",help="the file, must be headered as \"Pathway-id	Pathway-name	Gene-id	Gene-name\"",required=True)

args = parser.parse_args()
vcffi = args.fi

full = csv.DictReader(open(vcffi,'r'),delimiter="\t")


#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*

####def anyNone(rets):

def getGenes(pathid,pth):#idea here, get a gene by position, and step forward only.

    count = 0
    
    (pwyid,pwyname) = pathid.split(':')
    while count < len(pth):

        frontgene = pth[count]

        for genes in pth[count + 1:len(pth)]:

            (geneid,genename) = genes.split(":")
            (frontid,frontname) = frontgene.split(":")

            print(pwyid + "\t" + pwyname + "\t"  + geneid + "\t" + genename + "\t" + frontid + "\t" + frontname )

        count = count + 1



#---------------------------------main-----------------------------------#

pre_dict = {}


for line in full:#load dict ass array per pathway.
    pathkey = line['Pathway-id'] + ":" + line["Pathway-name"]
    if pathkey in pre_dict:

        if "unknown" not in line['Gene-id']:
            pre_dict[pathkey].append(line['Gene-id'] + ":" + line['Gene-name'])

    else:
        pre_dict[pathkey] = []


print('Pathway-id\tPathway-name\tGene-id\tGene-name\tTarget-id\tTarget-name')

for path in pre_dict:

    getGenes(path,pre_dict[path]) 
    










