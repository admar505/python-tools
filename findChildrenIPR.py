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


#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*
def anyNone(rets):

    size = len(rets)
    none_ct = size  #set to size, reduce as nones go


def getChildren(rootelement):
    print(rootelement.attrib['id'])
    for child in rootelement.findall('./channel/item'):
        print(child)
        #print(child.attrib['db_xref'])


    return None


#####--------------------------------MAIN-------------------------------------------######

for element in iprtree.getroot():
    try:

        if element.attrib['id'] in taglst:

            print(element.attrib['id'])
            getChildren(element)


    except KeyError:
        
        next

        #print("unknown root")





