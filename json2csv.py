#!/usr/bin/python
import json
import argparse

parser = argparse.ArgumentParser(description="convert PGX json to csv")
parser.add_argument("--file",help="json from PGX")
args = parser.parse_args()
#pull in file
infi = args.file

print infi
data_file = open(infi,"r")
data_json = json.load(data_file)
#print data_json['groups']
print data_json.keys()
#print json.dumps(data_json)
print "relatedChemicals:" + str(data_json['relatedChemicals'][0]['id'])
#chemical_id = data_json['id']
#for i in chemical_id:
#    print i


#print "haplotype\tassoc.strength\tphenotype\tannot\t"
for element in data_json['groups']:
    for i in element['genotypes']:
#        print element.keys()
#    for ele in element['annotations']:
#       print str(ele['markdown']['markdown']).replace('\n',' ')
        print str(data_json['relatedGenes'][0]['symbol']) + "\t" +  i + "\t" +   str(data_json['relatedChemicals'][0]['id']) +  "\t"  + str(element['strength']['term']) + "\t"  +  str(element['name']).strip()   + "\t" + str(element['annotations'][0]['markdown']['markdown']).replace('\n',' ') + "\t"  +   str(element['annotations'][1]['markdown']['markdown']).replace('\n',' ')  + "\t" + str(element['annotations'][2]['markdown']['markdown']) + "\t" + str(element['annotations'][3]['markdown']['markdown']).replace('\n',' ') + "\t" + str(element['annotations'][4]['markdown']['markdown']).replace('\n',' ')

        #for ann in 1:len(element['annotations']):
        #    print str(element['annotations'][ann])
#print data_json['variants'] #this does not make sense.
#for ele2 in data_json['markdown']:
    #print ele2["id"]
#    for eledeep in ele2:
#        print eledeep
   # for ele2_more in ele2.iterkeys():
    #    print ele2_more

