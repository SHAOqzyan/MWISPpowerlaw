#from madda import  myG210
from astropy.table import Table
from powerLawComponent1 import powerLaw1
import numpy as np
#doG210 = myG210()
import glob
doPowerLaw=powerLaw1()


#first compare the


class testPL:

    saveAlphaPath = "./savePowerlawPhysicalArea/"
    saveTBPath= "./saveTBPath/"
    con1G2650= np.arange( 4 , 8, 1)
    con2G2650= np.arange( 8 , 20, 1 ) #start from 5
    con3G2650= np.arange( 11 , 28, 1) #start from 6

    conTypeG2650={1:con1G2650,2:con2G2650,3:con3G2650}
    cutoffList= [2, 2.5,3.0,3.5,4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
    rmsCO12=0.5
    def __init__(self):
        pass

    def comparePwerLaw(self,verbose=True):
        """
        compare the result of two script
        :return:
        """
        minArea=0.03
        #testTBFile= "/home/qzyan/WORK/diskMWISP/fillingFactorData/tmpFiles/rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"

        testTBFile= "rawLocalCO12rawLocalCO12_SmFactor_1.0_NoiseAdd_0.0dbscanS2P4Con1_Clean.fit"

        tb=Table.read(testTBFile)


        areaArray = tb["area_exact"]/3600.

        maxArea=  np.max(areaArray)


        rawArea = areaArray[areaArray>=minArea] #



        print "???????????????????????????????????????????????????????"
        areaError=0.2*rawArea

        meanAlpha, stdAlpha = doPowerLaw.getAlphaWithMCMCWithErrorMultiChains(rawArea,areaError, minArea=minArea, maxArea=maxArea)
        #part1 = doPowerLaw.fitPowerLawWithMCMCcomponent1(rawArea,  minV=minArea, maxV=maxArea)


        print meanAlpha,stdAlpha

    def removeAllEdges(self, TBList):
        """

        :param TBList:
        :return:
        """
        newList = []

        for eachTB in TBList:
            newList.append(self.removeWrongEdges(eachTB))

        return newList
    def removeWrongEdges(self, TB):

        if TB == None:
            return None
        processTB = TB.copy()

        # remove cloudsThat touches the noise edge of the fits

        # part1= processTB[ np.logical_and( processTB["x_cen"]>=2815 ,processTB["y_cen"]>= 1003  )   ] #1003, 3.25

        # part2= processTB[ np.logical_and( processTB["x_cen"]<= 55 ,processTB["y_cen"]>= 1063  )   ] #1003, 3.25

        if "peak" in TB.colnames:  # for db scan table

            part1 = processTB[np.logical_or(processTB["x_cen"] > 26.25, processTB["y_cen"] < 3.25)]  # 1003, 3.25
            # part1= processTB[ np.logical_or( processTB["x_cen"]>26.24166667 ,processTB["y_cen"] < 3.25833333 )   ] #1003, 3.25

            part2 = part1[np.logical_or(part1["x_cen"] < 49.25, part1["y_cen"] < 3.75)]  # 1003, 3.25
            # part2= part1[ np.logical_or( part1["x_cen"]<49.24166667 ,part1["y_cen"]<  3.75833333 )   ] #1003, 3.25

            return part2
        else:  # dendrogram tb

            part1 = processTB[np.logical_or(processTB["x_cen"] < 2815, processTB["y_cen"] < 1003)]  # 1003, 3.25

            part2 = part1[np.logical_or(part1["x_cen"] > 55, part1["y_cen"] < 1063)]  # 1003, 3.25

            return part2

    ##########################################################################################
    def selectTBFormal(self,TB,cutOff=2 ,  pixN=16, minDelta= 3 ,hasBeam=True,minChannel=3 ,verbose=False,removeEdge=True ):
        """
        # This is the most strict critera to select, first by pixN, second by miNDelta, which only affect peak, the peak is calculated by cuOff+minDelta,
        # minChannel and has Beam would also be done
        :param TB:
        :param areaPix:
        :param conChannel:
        :return:
        """

        #first, check cuOff, to prevent wrong cutOff

        #first voxel
        filterTB = TB[TB["pixN"] >= pixN ]

        #second peak
        filterTB = filterTB[filterTB["peak"] >= (minDelta+cutOff)*self.rmsCO12]

        #third by beam,
        if hasBeam: # this is
            filterTB = filterTB[filterTB["has22"] >=  0.5  ]

        #select by minCHannel
        filterTB = filterTB[filterTB["allChannel"] >= minChannel  ]


        #reamove edged
        if removeEdge:
            tmpList=self.removeAllEdges([filterTB])

        else:
            return filterTB


        return tmpList[0]




    def getTBByCutOffList(self, cutOffList, minPts, conType, calCode="G2650Local" , getCleanLW=True ,selectionCode="selectFormal"  ):
        """

        :param cutOffList:
        :param minPts:
        :param conType:
        :param calCode:
        :return:
        """
        #
        returnTBList= []
        searchPath=self.saveTBPath #"/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/"
        for eachCutoff in cutOffList:

            searchStr= searchPath+"{}dbscanS{}P{}Con{}.fit".format( calCode,  eachCutoff  ,  minPts ,    conType  )
            if getCleanLW:
                searchStr = searchPath + "{}dbscanS{}P{}Con{}*CleanWithLW.fit".format(calCode, eachCutoff, minPts, conType)

            a= sorted( glob.glob(  searchStr  ) ) #sort this in names

            if len(a)==0:
                returnTBList.append(None)

            else:

                tmpTB = Table.read( a[0] )

                #tmpTB=self.selectTB(tmpTB,areaPix=areaPix,conChannel=conChannel) #by pixel
                #tmpTB=self.selectTBByPeak(tmpTB,minDelta=eachCutoff+3 ) #by peak
                tmpTB=self.selectTBFormal(tmpTB,cutOff= eachCutoff )

                returnTBList.append( tmpTB )


        #remove wrong edges
        returnTBList=self.removeAllEdges(returnTBList)

        return returnTBList

    def getPhysicalAlphaList(self, TBList):
        """

        :param TBList:
        :return:
        """

        #####

        alphaList = []
        alphaErrorList = []

        # physicalEdges = np.linspace(0, 100, 1000)  # square pc^2
        # physicalCenter = self.getEdgeCenter(physicalEdges)
        length = 1500 * np.deg2rad(0.5/60)
        compoleteAreaPhysical = length ** 2 * 4  # 4 pixels

        for eachTB in TBList:
            # realArea = doDBSCAN.getRealArea(eachTB)
            print "Total number of clouds, ", len(eachTB)
            phyiscalArea, physicalAreaError = self.getPhyscialAreaAndError( eachTB )

            # binN, binEdges = np.histogram(realArea, bins=physicalEdges)

            # calculate alpha

            # meanA, stdA = doDBSCAN.getAlphaWithMCMC(realArea, minArea=compoleteAreaPhysical, maxArea=None, physicalArea=True)
            meanA, stdA = doPowerLaw.getAlphaWithMCMCWithErrorMultiChains(phyiscalArea, physicalAreaError,    minArea=compoleteAreaPhysical,  maxArea=np.max(phyiscalArea))

            alphaList.append(meanA)

            alphaErrorList.append(stdA)

        return alphaList, alphaErrorList



    def getPhyscialAreaAndError(self,TB):
        """

        :param TB:
        :return:
        """

        # print TB.colnames

        processTB=TB.copy()
        v = processTB["v_cen"]
        dis = (0.033 * v + 0.180) * 1000  # pc
        processTB["fakeDistance"] =     dis

        # dis= ( 0.033*v + 0.175)*1000 # pc
        #

        selectionCriteria= np.logical_and(dis>0, dis<=1500   )

        goodClouds= processTB[selectionCriteria] ###

        NArray = goodClouds["area_exact"] / 0.25

        # print N,eachR["pixN"]
        length = goodClouds["fakeDistance"] * np.deg2rad(0.5 / 60)  # square pc
        physicalArea = length ** 2 * NArray # eachR["pixN"]  #*10000

        dError=161 #pc
        physicalAreaError=2* goodClouds["fakeDistance"]*dError*(np.deg2rad(0.5 / 60) )**2*NArray

        return physicalArea,physicalAreaError


    def calPhysicalAlpha(self,conType=1):
        """

        :return:
        """

        for eachminPts in self.conTypeG2650[conType]:

            tbList = self.getTBByCutOffList(self.cutoffList,eachminPts,conType)


            #######
            alphaArea, alphaAreaError = self.getPhysicalAlphaList(tbList)

            print "PhysicalAreaMinPts{}Con{}".format( eachminPts, conType)
            print alphaArea,alphaAreaError

            saveTagAlpha=self.saveAlphaPath+"PhysicalAreaMinPts{}Con{}".format( eachminPts, conType)
            saveTagAlphaError=self.saveAlphaPath+"PhysicalAreaErrorMinPts{}Con{}".format( eachminPts, conType)

            np.save(saveTagAlpha, alphaArea  )
            np.save(saveTagAlphaError, alphaAreaError  )


    def calMassAlpha(self,conType=1):
        """

        :return:
        """

        for eachminPts in self.conTypeG2650[conType]:

            tbList = self.getTBByCutOffList(self.cutoffList,eachminPts,conType)


            #######
            alphaMass , alphaMassError = self.getPhysicalAlphaList(tbList)

            print "massAlphaMinPts{}Con{}".format( minPts, conType)
            print alphaArea,alphaAreaError

            saveTagAlpha = self.saveAlphaPath+"massAlphaMinPts{}Con{}".format( eachminPts, conType)
            saveTagAlphaError = self.saveAlphaPath+"massAlphaMinPts{}Con{}".format( eachminPts, conType)

            np.save(saveTagAlpha, alphaMass  )
            np.save(saveTagAlphaError, alphaMassError  )







    def ZZZ(self):
        pass


if 1:

    doPL=testPL()
    doPL.calPhysicalAlpha( conType=1 )
    doPL.calPhysicalAlpha( conType=2 )
    doPL.calPhysicalAlpha( conType=3 )