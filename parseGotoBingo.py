#!/usr/bin/env python3

import argparse
import string
import sys
import subprocess
import re

parser = argparse.ArgumentParser(description='Creates a BinGO format for use with Bingo (and other tools). Will parse into separate files, ')
parser.add_argument("--ann", help="The go.obo I have programmed this before, and will now do so again ",required=True)
parser.add_argument("--spc", help="The name of the org. Doesnt matter, but use it. ")



##4 files:
##namespace: biological_process
##namespace: cellular_component
##namespace: external
##namespace: molecular_function

args = parser.parse_args()
species = args.spc

#open the out files

bioproc = open('biological_process.txt','w')
cellcomp = open('cellular_component.txt','w')
ext = open('external.txt','w')
molfunc = open('molecular_function.txt','w')


######################----------------defs----------------#############################

def Printer(go,dsc,ais,part,out):
   

    isit = [] 
    prted = []

    isit.append(ais.keys())
    prted.append(part.keys())
    
    if len(ais.keys()) == 0 and len(part.keys()) == 0: 
        out.write(str(go) + " = " + str(dsc) + "\n")
        
    elif len(ais.keys()) == 0 and len(part.keys()) >= 1:
        out.write(str(go) + " = " + str(dsc) + " [partof: " + str(' '.join(part.keys())) + " ]\n")


    elif len(ais.keys()) >= 1 and len(part.keys()) >= 1:
        out.write(str(go) + " = " + str(dsc) + " [isa: "  + str(' '.join(ais.keys()))    + " ] [partof:" + str(' '.join(part.keys())) + " ]\n")


    elif len(ais.keys()) >= 1 and len(part.keys()) == 0:
        out.write(str(go) + " = " + str(dsc) + " [isa:" + str(' '.join(ais.keys())) + " ]\n")


######################----------------main----------------#############################
#gohash = {} #storing the stuff.

obofi = open(args.ann,'r')

bioproc.write("(type=Biological Process)(curator=Invaio)\n")
cellcomp.write("(type=Biological Process)(curator=Invaio)\n")
ext.write("(type=Biological Process)(curator=Invaio)\n")
molfunc.write("(type=Biological Process)(curator=Invaio)\n")

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


