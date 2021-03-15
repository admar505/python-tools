#!/usr/bin/env python3
import argparse
import subprocess
import os
import re
import csv
parser = argparse.ArgumentParser(description='Given the NovoGene qc file, create the directory of veribench import info')
parser.add_argument("--batch", help="Batch name, just a string",required=True)
parser.add_argument("--date",help="date to use",required=True)
parser.add_argument("--qc", help="The QC file from the novogene report",required=True)

args = parser.parse_args()
date = args.date
batchname = args.batch
qc_handle = open(args.qc,'r')
qcfile = csv.DictReader(qc_handle,delimiter="\t")

##newer version, for version 2 of novogene qual metrics.

#<---userdefs--->##

def getflowcell(filedict):

    firstline = next(filedict)
    #flowcell = firstline["Lane"].split('_')[0]

    flowcell = firstline["Lane"]
    return flowcell

def qcwrite(out,qin):
    out.writeheader()

    def getlane(lanein):
        laneout = str(lanein).split('_')
        lanenum = str(laneout[1]).split('L')[1]
        return(lanenum)

    for q in qin:#each line, write to the file and stuff
        #lanes = getlane(q['Lane'])
        lanes = 1
        out.writerow({'Sample Name':q['Sample name'],'% >= Q30 bases':q['Q30(%)'], 'Lane':lanes})

def clustwrite(cout,qin):

    for q in qin:
        cout.write(q['Sample name'] +","+ q['Raw data'] +"\n")


def indexwrite(iout,qin):
    for q in qin:
        iout.write(q['Sample name'] +","+ q['Raw reads'] +"\n")


# ------main------

##get ids together.
os.makedirs(batchname,exist_ok=True)
os.makedirs(batchname + "/runInfo",exist_ok=True)
os.makedirs(batchname + "/qc",exist_ok=True)

clusterfilename = batchname + "/runInfo/" +  str(date) + "_NOVOGENERUN_" + getflowcell(qcfile) + ".cluster_metrics_summary.csv"
indexfilename = batchname  + "/runInfo/" + str(date) + "_NOVOGENERUN_" + getflowcell(qcfile) + ".index_metrics_summary.csv"
qcfn = batchname + "/qc/qc30_metrics.csv"
qc_handle.seek(1)

qcfileout = open(qcfn,'w')
clusterfile = open(clusterfilename, "w" )
indexfile = open(indexfilename, "w")

#make qc
qcfields = ['Sample Name', '% >= Q30 bases',  'Lane']
qcout = csv.DictWriter(qcfileout, fieldnames=qcfields,delimiter="\t")

qc_handle.seek(0)
qcwrite(qcout,csv.DictReader(qc_handle,delimiter='\t'))

qc_handle.seek(0)
clustwrite(clusterfile,csv.DictReader(qc_handle,delimiter='\t'))

qc_handle.seek(0)
indexwrite(indexfile,csv.DictReader(qc_handle,delimiter='\t'))






