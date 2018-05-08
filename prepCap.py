#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import vcf
from subprocess import call
import csv
import random
sys.path.append("/home/nucleo/lib/PyVGRes")
import vgr

parser = argparse.ArgumentParser(description="select infor from VCF")
parser.add_argument("--vcf",help="vcf file, bgzipped and tabixed as -p vcf",required=True)
parser.add_argument("--list",help="lost to search through, in CAP format",required=True)

args = parser.parse_args()
vcffi = args.vcf
varslst = args.list

vcf_full = vcf.Reader(open(vcffi,'r'))
necessary_vars = csv.DictReader(open(varslst,'r'),delimiter="\t")

#parse results in a map or dict, or what??

#-------------------------------------here by DEFSgONS!!----------------------------------*
def recallVar(fudge):
    print "FUCK"



def writeFinal():
    return None



def mergeNWrite(vgrfile,caps):
    complete = vgr.Reader(open(vgrfile,'r'))





    return None


def kickAnnotations(vcfln,capln,vcf_4_header):
    returnvals = []
    #command = "/vbin/GoPipeRUN/goVarAnnotateAndID.sh" + rsid
    #call(["/vbin/GoPipeRUN/goVarAnnotateAndID.sh",rsid])
    rsid = capln['Chr'] + "."  + capln['Pos']
    vcf_write = vcf.Writer(open(rsid + '.Merged.vcf','w'),vcf_4_header)
    vcf_write.write_record(vcfln)
    sys.stdout.flush()

    effout = open(rsid + '.efd.vcf',"w")
    errout = open(rsid + '.stderrr',"w")

    call(['java -jar /vbin/snpEff/snpEff.jar eff -i vcf -csvStats -hgvs hg19 ' + rsid + '.Merged.vcf'],shell=True,stdout=effout,stderr=errout)
    sys.stdout.flush()
    effout.close()

    call(['/vbin/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl --force_overwrite --vcf -i ' + rsid + '.Merged.vcf -o ' + rsid + '.VEP.vcf --dir /vbin/ensembl-tools-release-86/  --everything --species homo_sapiens --cache --refseq --offline --assembly GRCh37 --fasta /vbin/ref/hg19.nochr.fa'],shell=True,stderr=errout)
    sys.stdout.flush()
    complete = open(rsid + '.COMPLETE.txt',"w")

    call(['/vbin/Perl/matchPathsAndMergeCallers.2.pl  --evs /ref/EVS_AF.vcf -p /ref/FULL.GENOME.lookUP.txt -eff ' + rsid + '.efd.vcf -vep ' + rsid + '.VEP.vcf -trn /ref/WG_good_genes.lst  -exac /ref/ExAC.r0.3.1.sites.af.vcf  -hgmd /ref/Homo_sapiens.HGMD.hg19.chr.vcf'],shell=True,stdout=complete,stderr=errout)

    mergeNWrite(rsid + '.COMPLETE.txt',capln)






def getVCFLine(varln,vfile):

    pos = re.search('(\d+)\-(\d+)',varln['Pos'])

    if pos is not None:
        found_vcf = vfile.fetch(varln['Chr'],int(pos.group(1)) - 1,int(pos.group(2)) + 1)

    else:
        found_vcf = vfile.fetch(varln['Chr'],int(varln['Pos']) -1,int(varln['Pos']) + 1)

    return_line = None

    for captured in found_vcf:
        return_line = captured

    return return_line




#####----------------MAIN--------------####      #####----------------MAIN--------------####



for needed in necessary_vars:
    extracted = getVCFLine(needed,vcf_full)
    if extracted is not None:
        kickAnnotations(extracted,needed,vcf_full)


    else:
        recallVar(needed)





















