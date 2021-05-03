#!/usr/bin/python
import sys,os,re,fileinput,argparse
import vcf
from subprocess import call
#idea: grab the full vcf in a dict/map; grab the omicia vcf in a dict/map; 
#then, then make a map to of res. read through the vcf, and add to reslines. If res in vcf, then 
#add score and rsid to EFF_EFFECT. if not, then grab line from big vcf, and process like olden times as if new. call recover.RESULTS.txt
parser = argparse.ArgumentParser(description="Make the file unique, if there are multiple lines that are different.")
parser.add_argument("--fi",help="The, uh, file.")


args = parser.parse_args()
fifi = args.fi
#special load for vcffi? (vcffi,"r")

keyid = {}

header = ""


with open(fifi) as in_file:

    for line in in_file:
        if "category" in line:
            header = line.rstrip() 

        lns = line.split("\t")
        
        if str(lns[3]) in keyid:
            idkey = keyid[lns[3]].split("\t")
            
            if int(idkey[6]) <= int(lns[6]):
                keyid[lns[3]] = line.rstrip()

        else:
            keyid[lns[3]] = line.rstrip()
        

print(header)

for prln in sorted(keyid.keys(),reverse=True):

    print(keyid[prln])
    
    


