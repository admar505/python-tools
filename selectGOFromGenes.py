#!/usr/bin/env python3

import argparse
import string
import sys
import subprocess
import re

parser = argparse.ArgumentParser(description='Takes in a list of genes with associated GO assignments. This is related to the script parseGotoBingo.py.')
parser.add_argument("--ann", help="The go.obo I have programmed this before, and will now do so again ",required=True)
parser.add_argument("--spc", help="GO annotated gene ids. so, like this, you see NP_001153504.1 has   0005887, that a tab in there you see.")



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

## Part way there.
##But NEED to get info on all ids, where to get that? probably from the go.obo, and I need to 
##get all descs in. so, some better way of...
##well, ok, NAh, I need to get this into one sing file, and descs and all. its in the obo file.
##
##Parse file twice, with the second parse is the stuff I already have. 
##first pass:
##  1. get the id to desc.
##  2. get the pairwise map, recurse on that.
##
##
######################----------------defs----------------#############################


def recurseGraph(gid,gdesc,grmap,out,dontgoagain):
    
    if str(gid) in grmap and str(gid) not in dontgoagain:
        
        out.write(gid + "\t" + str(gdesc[gid]) + "\tdescendant\t" + str(grmap[gid]) + "\t" + str(gdesc[grmap[gid]]) + "\n"  )

        dontgoagain[str(gid)] = gid
        dontgoagain[str(grmap[gid])] = grmap[gid]

        recurseGraph(grmap[gid],gdesc,grmap,out,dontgoagain)


    


def Printer(go,dsc,ais,part,out,gdsc,grpaths,been):##NEW, recurse and pull the descs on print.
   
        
    if len(ais.keys()) == 0 and len(part.keys()) >= 1:

        for gothing in part.keys():
            out.write(str(go) + "\t" + str(dsc)  +  "\tpartof\t" + str(gothing) + "\t" + str(gdsc[gothing]) + "\n")
            recurseGraph(go,gdsc,grpaths,out,been)

    elif len(ais.keys()) >= 1 and len(part.keys()) >= 1:

        for gothing in ais.keys():
            out.write(str(go) + "\t" + str(dsc) +  "\tisa\t"  + str(gothing) + "\t" + str(gdsc[gothing]) +  "\n")

            recurseGraph(go,gdsc,grpaths,out,been)

        for gopart in part.keys():
            out.write(str(go) + "\t" + str(dsc) + "\tpartof\t"  + str(gopart) +  "\t" + str(gdsc[gopart]) + "\n")

            recurseGraph(go,gdsc,grpaths,out,been)

    elif len(ais.keys()) >= 1 and len(part.keys()) == 0:

        for gothing in ais.keys():
            out.write(str(go) +  "\t" +  str(dsc) + "\tisa\t" + str(gothing) + "\t"  +  str(gdsc[gothing]) + "\n")
            
            recurseGraph(go,gdsc,grpaths,out,been)
            
            
            
            

######################----------------main----------------#############################

###  FILE handling ###

#gohash = {} #storing the stuff.

firstpassfi = open(args.ann,'r')
graph = {}##the tree structure

auxgodesc = {}##just other descriptions that belong to the go terms
seen = {}#this is term condition on the recurse.

base_go_id = ""

for fpl in firstpassfi:


    if "id: GO" in fpl and "alt_id:" not in fpl:
        #print(line)
        goidre = re.search("GO\:(\d+)",fpl)
        base_go_id = goidre.group(1)


    elif "name" in fpl:
        graph[base_go_id] = fpl.split(":")[1].strip()
        auxgodesc[base_go_id] =  fpl.split(":")[1].strip()


    elif "is_a: GO" in fpl or "intersection_of:" in fpl:
        rs = re.search("GO\:(\d+)",fpl)
        rsid = rs.group(1)

        graph[rsid] = base_go_id
        graph[base_go_id] = rsid

        dscr = re.search("! (\S+.*)",fpl)
        auxgodesc[rsid] = dscr.group(1)
        


    elif "part_of GO" in fpl:
        rs = re.search("GO:(\d+)", fpl)

        graph[base_go_id] = rs.group(1)
        graph[rs.group(1)] = base_go_id
        
        dscr = re.search("! (\S+.*)",fpl)
        auxgodesc[rsid] = dscr.group(1)




firstpassfi.close()


obofi = open(args.ann,'r')
saveGO = {}## dictionary of GO ids to save.


genefi = open(args.spc,'r')

bioproc.write("id\tid_desc\trelation\ttarget\ttarget_desc\n")
cellcomp.write("id\tid_desc\trelation\ttarget\ttarget_desc\n")
ext.write("id\tid_desc\trelation\ttarget\ttarget_desc\n")
molfunc.write("id\tid_desc\trelation\ttarget\ttarget_desc\n")

for gln in genefi:
    (geneid,haser,GOid) = gln.split("\t")
    saveGO[GOid.strip()] = geneid.strip()

    
    bioproc.write(geneid + "\tgene_of_unknown_function\thasa\t" + str(GOid.strip()) + "\t" + str(auxgodesc[GOid.strip()]) + "\n" )
    cellcomp.write(geneid + "\tgene_of_unknown_function\thasa\t" + str(GOid.strip()) + "\t" + str(auxgodesc[GOid.strip()]) + "\n")
    ext.write(geneid + "\tgene_of_unknown_function\thasa\t" + str(GOid.strip()) + "\t" + str(auxgodesc[GOid.strip()]) + "\n")
    molfunc.write(geneid + "\tgene_of_unknown_function\thasa\t" + str(GOid.strip()) + "\t" + str(auxgodesc[GOid.strip()]) + "\n")

    
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
                Printer(goid,desc,isa,partof,bioproc,auxgodesc,graph,seen)

            elif "cellular_component" in namespace:
                Printer(goid,desc,isa,partof,cellcomp,auxgodesc,graph,seen)
        

            elif "external" in namespace:
                Printer(goid,desc,isa,partof,ext,auxgodesc,graph,seen)

            elif "molecular_function" in namespace:
                Printer(goid,desc,isa,partof,molfunc,auxgodesc,graph,seen)
 


        goid = ""
        desc = ""
        namespace = ""


