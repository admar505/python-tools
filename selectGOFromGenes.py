#!/usr/bin/env python3

import argparse
import string
import sys
import subprocess
import re

parser = argparse.ArgumentParser(description='Takes in a list of genes with associated GO assignments. This is related to the script parseGotoBingo.py.')
parser.add_argument("--ann", help="The go.obo I have programmed this before, and will now do so again ",required=True)
parser.add_argument("--spc", help="GO annotated gene ids. so, like this, you see NP_001153504.1    0005887, that a tab in there you see.")


args = parser.parse_args()
species = args.spc

#open the out files

bioproc = open('biological_process.txt','w')
cellcomp = open('cellular_component.txt','w')
ext = open('external.txt','w')
molfunc = open('molecular_function.txt','w')



##So this is new, and the point is to trace/augment the output file with the 
##etra nodes that are the is_a or part_of as the join.
##make that the interaction type. TRY IT!!!
##use the spc annotated list as the key, if the id is in, then save all the stuff, 
##I think it is a quick wrapper call for is it?? in??


######################----------------defs----------------#############################

def Printer(go,dsc,ais,part,out):
   

    isit = [] 
    prted = []

    isit.append(ais.keys())
    prted.append(part.keys())
    
        
    if len(ais.keys()) == 0 and len(part.keys()) >= 1:

        for gothing in part.keys():
            out.write(str(go) + " partof " + str(gothing) + "\n")


    elif len(ais.keys()) >= 1 and len(part.keys()) >= 1:

        for gothing in ais.keys():
            out.write(str(go) + " isa "  + gothing + "\n")

        for gopart in part.keys():
            out.write(str(go) + " partof "  + gothing + "\n")

    elif len(ais.keys()) >= 1 and len(part.keys()) == 0:

        for gothing in ais.keys():
            out.write(str(go) + " isa " + gothing + " ]\n")


######################----------------main----------------#############################

###  FILE handling ###

#gohash = {} #storing the stuff.

obofi = open(args.ann,'r')

bioproc.write("(type=Biological Process)(curator=Invaio)\n")
cellcomp.write("(type=Biological Process)(curator=Invaio)\n")
ext.write("(type=Biological Process)(curator=Invaio)\n")
molfunc.write("(type=Biological Process)(curator=Invaio)\n")

saveGO = {}## dictionary of GO ids to save.


genefi = open(args.spc,'r')

for gln in genefi:
    (geneid,GOid) = gln.split("\t")
    saveGO[GOid.strip()] = geneid.strip()

    
genefi.close()


isa = {}
partof = {}
goid = "" 
desc = ""
namespace = ""

for line in obofi:
    if "[Term]" in line:    #start collecting
                
        isa = {}#           the isa term holder.
        partof = {}#        hold part of, put after the isa ALWAYS

    elif "id: GO" in line and "alt_id:" not in line:
        #print(line)
        goidre = re.search("GO\:(\d+)",line)
        goid = goidre.group(1)

    elif "namespace" in line:
        namespace = line.split(":")[1].strip()


    elif "name" in line:
        desc = line.split(":")[1].strip()
        
    elif "is_a: GO" in line:
        rs = re.search("GO\:(\d+)",line)
        isa[rs.group(1)] = "isa"

    elif "part_of GO" in line:
        rs = re.search("GO:(\d+)",line)
        partof[rs.group(1)] = "part_of"

    elif len(line.strip()) == 0:#sort into the file type
        if goid in saveGO:#this checks to see if the base go is in here. problem, I dont know
                          #if this will save higher ids.

            if "biological_process" in namespace:
                Printer(goid,desc,isa,partof,bioproc)

            elif "cellular_component" in namespace:
                Printer(goid,desc,isa,partof,cellcomp)
        

            elif "external" in namespace:
                Printer(goid,desc,isa,partof,ext)

            elif "molecular_function" in namespace:
                Printer(goid,desc,isa,partof,molfunc)
 


        goid = ""
        desc = ""
        namespace = ""


