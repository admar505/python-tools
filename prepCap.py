#!/usr/bin/env python
import sys,os,re,fileinput,argparse
import vcf
import csv
import random
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
def anyNone(rets):

    size = len(rets)
    none_ct = size  #set to size, reduce as nones go


    for tagkey in rets:
        checkfornone = re.match('.*?\[None\].*',str(rets[tagkey]))#LOGIC: if there is a none check for none will be NOT NONE??

        try:

            if checkfornone is  None:
                none_ct = none_ct - 1 #decided to use as --, as it is more sensible.

        except (AttributeError, IndexError) as e:
            dn = open(os.devnull,'w')

    final_return = None

    if none_ct < size:##indicates that no NONE vals were found in entirety.
        final_return = rets

    return final_return


def getTags(tags, varset):
    retval = {}

    for tagval in tags:

        if tagval in varset:
            #retval.append(tagval +  '=' + str(varset[tagval]))
            retval[tagval] = varset[tagval]

    #make a loop or def() that checks of at least one is not none.

    return_final = anyNone(retval)
    #print return_final
    return return_final


def kickAnnotations(vcfln,capln):
    returnvals = []
    #command = "/vbin/GoPipeRUN/goVarAnnotateAndID.sh" + rsid
    #call(["/vbin/GoPipeRUN/goVarAnnotateAndID.sh",rsid])
    effout = open(rsid + '.efd.vcf',"w")
    call(['java -jar /vbin/snpEff/snpEff.jar eff -i vcf -csvStats -hgvs hg19 ' + rsid + '.Merged.vcf'],shell=True,stdout=effout)

    call(['/vbin/ensembl-tools-release-78/scripts/variant_effect_predictor/variant_effect_predictor.pl --force_overwrite --vcf -i ' + rsid + '.Merged.vcf -o ' +rsid + '.VEP.vcf --everything --species homo_sapiens --cache --refseq --offline --assembly GRCh37 --fasta /vbin/ref/hg19.nochr.fa'],shell=True)
    complete = open(rsid + '.COMPLETE.txt',"w")
    call(['/vbin/Perl/matchPathsAndMergeCallers.2.pl  --evs /ref/EVS_AF.vcf -p /ref/FULL.GENOME.lookUP.txt -eff ' + rsid + '.efd.vcf -vep ' + rsid + '.VEP.vcf -trn /ref/WG_good_genes.lst  -exac /ref/ExAC.r0.3.1.sites.af.vcf  -hgmd /ref/Homo_sapiens.HGMD.hg19.chr.vcf'],shell=True,stdout=complete)





    return "\t".join(returnvals)



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
        kickAnnotations(extracted,needed)


    else:
        recallVar(needed)





















