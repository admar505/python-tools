#!/usr/bin/python
import sys,os,re,fileinput,argparse
import csv

parser = argparse.ArgumentParser(description="find var lines in answer bed that are for genes on Crick")
parser.add_argument("--tacos",help='answer bed to search through',required=True)

args=parser.parse_args()


def getHGVS(hgvs):
    get = re.search('.*?[0-9]([^0-9])\>(\w+)',hgvs)
    return(get.group(1) + get.group(2))



def getRSID(rs):
    ref = rs.split("-")[1]
    alt = rs.split("-")[2]

    return(ref + alt)

csvfi = csv.DictReader(open(args.tacos,'r'),delimiter='\t')

for var in csvfi: #parse the rsid and compare with hgvs call

    rsid = getRSID(var['rsid'])
    homogvs = getHGVS(var['homohgvs'])


    #print(rsid + "\t"  +homogvs)

    if rsid != homogvs:
        print(var['rsid'] +"\t"+ var['homohgvs'] + "\t" + var['hethgvs'] +"\t"+ var['homourl'] +"\t"+ var['heturl'] +"\t"+ var['wturl'])


