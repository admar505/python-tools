import sys,re,operator
import vcf
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

    def __init__(self,variants,rptfile=None):#I need to load the above to here, how doI do that.
        #print "alive!---------------------------------------"
        alts = loadAlts()                   #loading...
        self.alts = alts.altLoader(variants)#all the alts
        self.WT = variants.REF              #ref, for fast pulling
        self.GL = variants.samples[0]['GL'] #likelihoods.
        self.QUAL = variants.QUAL           #qual
        self.full = variants                #the whole thing for everything else
        self.rpt = rptfile

    @property
    def test(self):
        return "alive"

    @property
    def retQUAL(self):
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
        #index += 1
        return self.alts[index]

    def __retWTorALT__(self,posA,posB):#returns GT for numerical position input
        thisGT = self.retREF +  ","  + self.retREF#create dictionary that both can pull from

        if posA == 0 and posB > 0:
            thisGT = str(self.retREF) + "," + str(self.__retGT__(posB))

        elif posA > 0 and posB > 0:
            thisGT = str(self.__retGT__(posA)) + ","  + str(self.__retGT__(posB))

        return thisGT




    def __retGL__(self,numpos): #will return the GL given a pos
        GL_wanted = self.GL     #ok, hard here, sometimes array, sometimes not. so,
                                #have to adjust how to acess.
        if type(self.GL) != type(float()):#enter if not single val
            GL_wanted = self.GL[numpos]
            #for gl in self.GL:
                #GL_wanted = gl
        return GL_wanted


    @property
    def assGT_GL(self):     #to connect the GL and call pair.
                            #formula:F(j/k) = (k*(k+1)/2)+j from vcf4.2 doc.
                            #for 00,01,11,02,12,22
        gtglpairs = {}      ##$GT => GL is format.
        nuclA = 0
        nuclB = 0

        while nuclA + nuclB < len(self.alts) + len(self.alts):   #I think for this I want to do some crazy logic to assign positions.
            #nuclB = 0                                            #try:for increment B, then when A=B, set A to zero and increment
            checkpos  = str(nuclA) + str(nuclB)
            currentGT = self.__retWTorALT__(nuclA,nuclB)
            #GLpos = (nuclA * ((nuclA + 1)/2)) + nuclB
            GLpos = (nuclB * (nuclB + 1)/2) + nuclA
            currentGL = self.__retGL__(int(GLpos))
            #print str(currentGL) + "   POS IN GL:" + str(GLpos) + "\tcurrentGT:" + str(currentGT) + "\tposition call and check\t" + str(checkpos)
            gtglpairs[currentGT] = currentGL

            if nuclB == nuclA and nuclB + nuclA != 0:#this adjusts the GL timing.
                nuclA = 0
                nuclB += 1
            elif nuclA + nuclB == 0:
                nuclB += 1
            else:
                nuclA += 1

        return  gtglpairs

class detRepeats(object):

    def __init__(self,var_fi):#I need to load the above to here, how doI do that.
        #print "alive!---------------------------------------"
        #alts = loadAlts()                   #loading...
        self.full = var_fi                #the whole thing for everything else

    @property
    def test(self):
        return "alive"

    def __createHGVS__(self,unt,allA,allB,cdot):
        print "fudge"
        hgvs = "c."+ cdot + unt + "[" + allA + "]/[" + allB + "]"
        print hgvs
        return hgvs


    def __getRPT_GT__(self,unt,allA,allB,rptstart):
        alleleAcnt =  str(unt) * int(allA)# expand repeat
        allelA = rptstart + alleleAcnt
        alleleBcnt =  str(unt) * int(allB)
        allelB = rptstart + alleleBcnt
        hgvsgt = allelA + ',' + allelB
        if allA < allB:
            hgvsgt = allelA + ',' + allelB

        return hgvsgt

    def __cntLength__(self,vcfln,base,unit,wt,wtlen,cdot,rptstart):
        gt_and_hgvs = []

        if vcfln == 'WT+':
            gt_and_hgvs.append(self.__getRPT_GT__(unit,wtlen,wtlen,rptstart))
            gt_and_hgvs.append(self.__createHGVS__(unit,wtlen,wtlen,cdot))

        else:
            for alt in vcfln.ALT:
                print vcfln.samples[0]['GT']
                print base
                correct_gt_here = detGenoType(vcfln)    #do I need to test for het for sure? kinda.
                print correct_gt_here.assGT_GL#this gives me the GLs to sort through to get the GT if not WT+
                GTGL =  sorted(correct_gt_here.assGT_GL.items(), key=operator.itemgetter(1),reverse = True)[0]
                print GTGL
                alleles = str(GTGL)
                alleleA = alleles.split()[0]
                alleleB = alleles.split()[1]
                just_repeatA = re.search(r'%s(.*)' % base,str(alleleA))#need to find two lengths, alt and repeat.
                just_repeatB = re.search(r'%s(.*)' % base,str(alleleB))#need to find two lengths, alt and repeat.
                print str(alleleA) + "  what A"
                print str(alleleB) + "  what B"
                print just_repeatA.group(1)

                gt_and_hgvs.append(self.__getRPT_GT__(unit,wtlen,wtlen,rptstart))
                gt_and_hgvs.append(self.__createHGVS__(unit,wtlen,wtlen,cdot))

        return gt_and_hgvs



    def countRepeats(self,repeat):#give the repeat file.#example:  chr   start     stop     [minimum_base_unit]    [unit]  [wt+]
        rpt_size = None
        rvals = []

        for proto in repeat:
            rvals = proto.split()
            chrom = rvals[0]
            start = rvals[1]
            stop  = rvals[2]
            base  = rvals[3]
            rptfront = rvals[4]
            unit  = rvals[5]
            wt    = rvals[6]
            wtlen = rvals[7]
            cdot  = rvals[8]

        vcfsection = self.full.fetch(str(chrom),int(start),int(stop))
        #print vcfsection.call_rate()
        try:#(run info, lets do this)
            rpt_size = self.__cntLength__(next(vcfsection),base,unit,wt,wtlen,cdot,rptfront)

        except StopIteration:#return WT+ alleles.
            rpt_size = self.__cntLength__("WT+",base,unit,wt,wtlen,cdot,rptfront)


        return rpt_size

        #return self.rpt




#    @property
#    def returnAD(self):     #returns a dictionary of the allele ==> depth
#        alleleDepths = {}   #this is off a bit, I get error when I run and it doesnt make total sense,
#                            #weird contigency, if there are two alts and three AD vals, need to push the
#                            #second forward to take the third value, not the second. seems to be something to do with
#                            #the number '2' in that space and coupled (ie, adjacent) variants
#                            #potential bug here. delist this for now. check again with version 1.2
#        def __chooseD__(adindex,samples):
#            depthval = 0
#            try:#try catch here is to deal with the fact that AD is sometimes list
#                #and sometimes scalar val
#                depthval = samples[0]['AD'][adindex]
#            except(TypeError):
#                depthval = samples[0]['AD']
#            except(IndexError):
#                print "ERROR:AD index doesnt make complete sense in this context, it ran off end of AD" + str(adindex)
#            return depthval
#        ##------------
#        alleleDepths[str(self.WT)] = __chooseD__(0,self.full.samples)
#        try:#see note in __chooseD__ def, same problem here.
#            if len(self.alts.keys()) == 1 and len(self.full.samples[0]['AD']) == 3:#catching 1 alt and 3 AD contigency
#                print "ALT-INFO-ORIG:\t" + str(1) + "\t" + str(self.full.samples[0]['AD'])
#                alleleDepths[str(self.alts[1])] = __chooseD__(2,self.full.samples)
#            #    for alt in self.alts.keys():
#            #        print "ALT-INFO-EXC:\t" + str(alt) + "\t" + str(self.full.samples[0]['AD'])
#            #        alleleDepths[str(self.alts[alt])] = __chooseD__(alt,self.full.samples)
#        except:
#            #if len(self.alts.keys()) == 1 and len(self.full.samples[0]['AD']) == 3:#catching 1 alt and 3 AD contigency
#            #    print "ALT-INFO-ORIG:\t" + str(alt) + "\t" + str(self.full.samples[0]['AD'])
#            #    alleleDepths[str(self.alts[1])] = __chooseD__(2,self.full.samples)
#            #else:
#            for alt in self.alts.keys():
#                print "ALT-INFO-ELSE:\t" + str(alt) + "\t" + str(self.full.samples[0]['AD'])
#                alleleDepths[str(self.alts[alt])] = __chooseD__(alt,self.full.samples)
#        return alleleDepths




