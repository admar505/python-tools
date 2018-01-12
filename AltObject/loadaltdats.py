import sys,re
sys.path.append('./')
import altobject

class loadAlts(object):
    def __init__(self):
            self.alive = "I am alive!!"

    def altLoader(self,samples):#this loads the altobject.py
                                #strategy, give the number for each as it goes through. and it will return that var. ok. got it.
        alts = {}#dictionary
        genotype_index = 1

        for alt_choice in samples.ALT:
            #print alt_choice
            # myalts = altobject.AltObj(samples,genotype_index)

            alts[genotype_index] = alt_choice

            genotype_index += 1

        return alts

#here I am thinking sort of a subclass on this.
#currentalts = loadAlts(samples)
class detGenoType(object):

    def __init__(self,variants):#I need to load the above to here, how doI do that.
        print "alive!---------------------------------------"
        alts = loadAlts()
        self.alts = alts.altLoader(variants)
        self.WT = variants.REF
        self.GL = variants.samples[0]['GL']
        self.QUAL = variants.QUAL

    @property
    def test(self):
        return "alive"

    @property
    def retQUAL():
        return self.QUAL

    @property
    def retREF(self):
        return self.WT

    @property
    def defGT(self): #returns the sorted genotype.

        retGT = []
        for alt in self.alts:
            if self.alts[alt].amivalid == True:
                retGT.append(self.alts[alt].getcall)

        retGT.sort()

        return len(retGT)

    def __retGT__(self,index):#should return the alt at pos,
        index += 1
        return self.alts[index]

    def __retWTorALT__(self,posA,posB):#returns GT for numerical position input
        thisGT = self.retREF + self.retREF#create dictionary that both can pull from
        if posA == 0 and posB > 0:
            thisGT = str(self.retREF) + str(self.__retGT__(posB))
        elif posA > 0 and posB > 0:

            thisGT = str(self.__retGT__(posA)) + str(self.__retGT__(posB))

        return thisGT

    def __retGL__(self,numpos):#will return the GL given a pos
        pos = 0
        while pos < numpos:
            is_GL_wanted = self.GL[pos]#iterates through and finds the correct GL
            if pos == numpos:
                return is_GL_wanted
            pos += 1

    @property
    def assGT_GL(self):     #to connect the GL and call pair.
                            #formula:F(j/k) = (k*(k+1)/2)+j from vcf4.2 doc.
                            #for 00,01,11,02,12,22
        gtglpairs = {}      ##$GT => GL
        nuclA = 0
        #def __retCorrectGTforPos__()
        print self.GL
        while nuclA < len(self.alts):
            nuclB = 0

            while nuclB < len(self.alts):
                checkpos  = str(nuclA) + str(nuclB)
                print checkpos
                currentGT = self.__retWTorALT__(nuclA,nuclB)
                GLpos = nuclB * ((nuclB + 1)/2) + nuclA
                currentGL = self.__retGL__(int(GLpos))
                print str(currentGL) + "   POS IN GL:" + str(GLpos) + "\tcurrentGT:" + str(currentGT)
                gtglpairs[currentGT] = currentGL

                nuclB += 1
            nuclA += 1

        return  gtglpairs


