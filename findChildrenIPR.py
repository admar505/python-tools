#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import csv
import gzip
import random
import xml.etree.ElementTree as ET
parser = argparse.ArgumentParser(description="select member db ids for a given IPR id.")
parser.add_argument("--ipr",help="the IPR xml file",required=True)
parser.add_argument("--tag",help="tag to pull, add more --tag to pull more",required=True,action='append')

args = parser.parse_args()
iprfi = args.ipr
taglst = args.tag

ipr_full = gzip.open(iprfi,'r')
iprtree = ET.parse(ipr_full)

##extra iprs I wanna re-search through.send through a second time.
newtags = []


#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*

def iprCheck(doesthishaveipr):

    if "INTERPRO" in doesthishaveipr.attrib['db']:
        #print(str(doesthishaveipr) + "DEF::IPRCHECK::INSIDE")
        newtags.append(str(doesthishaveipr.attrib['dbkey']))


    
def getTags(parentIPR,dbxrefele):#$name of IPR (to print),  the dbxrefallele 

    dbs = ['TIGRFAMS','BLOCKS','CATH','CAZY','EC','CATHGENE3D','KEGG','PANDIT','PANTHER','PFAM','PRINTS','PRODOM','PROSITE','SCOP','TIGRFAMs','TIGRFAMS','SMART','INTERPRO']

    if dbxrefele.attrib['db'] in dbs:
        print(parentIPR  + "\t" + dbxrefele.attrib['db'] + "\t"  + dbxrefele.attrib['dbkey'])

        iprCheck(dbxrefele)




def childPrint(nameofIPR,node):

    for subnode in node.iter("db_xref"):
        getTags(nameofIPR,subnode)


def getChildren(iprname,rootelement):

    for child in rootelement:
        childPrint(iprname,child)        
  

def getParent(parentelem,iprs):

    try:
        if parentelem.attrib['id'] in iprs:
            #print(str(iprs) + "     the crazy stringer printer of iprs, did it add?")
            getChildren(parentelem.attrib['id'],parentelem)


    except KeyError:
        next


#####--------------------------------MAIN-------------------------------------------######

for element in iprtree.getroot():
    getParent(element,taglst)


#redo thesearch through the new collected IPR.

for additional in iprtree.getroot():
    getParent(additional,newtags)



