import sys,re
sys.path.append('./')
import altobject

class loadAlts(object):
    def __init__(self):
            self.alive = "I am alive!!"

    def altLoader(self,samples):#this loads the altobject.py
                                #strategy, give the number for each as it goes through. and it will return that var. ok. got it.
        alts = {}
        genotype_index = 1

            #print myalts.test()

        for alt_choice in samples.ALT:
            #print alt_choice
            myalts = altobject.AltObj(samples,genotype_index)
            alts[genotype_index] = myalts

            genotype_index += 1


            return alts

#here I am thinking sort of a subclass on this.
#currentalts = loadAlts(samples)
class detGenoType(object):

    def __init__(self,variants):#I need to load the above to here, how doI do that.
        #print "alive!"
        alts = loadAlts()
        self.alts = alts.altLoader(variants)
        self.WT = variants.REF
        self.GL = variants.samples[0]['GL']
    @property
    def test(self):
        return "alive"






    @property
    def defGT(self): #returns the sorted genotype.

        retGT = []
        for alt in self.alts:
            print alt
            if self.alts[alt].amivalid == True:
                retGT.append(self.alts[alt].getcall)



        retGT.sort()

        return len(retGT)


    @property
    def assocGTandGL(self): #formula:F(j/k) = (k*(k+1)/2)+j from vcf4.2 doc.
                            #for 00,01,11,02,12,22

        gtglpairs = {}
        nucl = 0
        #altypes = [self.WT,",".join(self.alts)]#add in the REF at pos zero.
        #print altypes
        test = ','.join(str(self.alts))
        for i in self.alts:
            print i
        while nucl < len(self.alts):
            vcfpos = nucl - 1

            if nucl == 1:

                dict_index =  altypes[vcfpos] + altypes[vcfpos - 1]

                gtglpairs[index] = self.GL[index]

            nucl += 1


        return self.GL




    @property
    def retREF(self):
        return self.WT
