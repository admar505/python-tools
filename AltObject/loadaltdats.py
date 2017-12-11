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

class determinGenoType(object):

        def __init__(self,alt):
            print "loading..."

