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
    def retREF(self):
        return self.WT
