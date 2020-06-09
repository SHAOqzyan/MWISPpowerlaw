
from progressbar import *
import numpy as np
from scipy.integrate import quad

from scipy.integrate import cumtrapz
import multiprocessing

def integrand(x, alpha,y,yError):
    """
    Does not need the normalize factor here, because the normal factor would bcalculated
    :param x:
    :param alpha:
    :param y:
    :param yError:
    :return:
    """
    #print x,alpha,y,yError

    #print x**(-alpha)*np.exp( - (y-x)**2*0.5/yError/yError )



    return x**(-alpha)*np.exp( - (y-x)**2*0.5/yError/yError )



class powerLaw1:
    thinning=10
    chainSample =  50
    burn_in=10
    nChain=40




    def __init__(self):
        pass



    # fitting function
    def calPowerLawProbLogcomponent1(self ,theta, dataArraym ,minV, maxV):
        alpha1 = theta[0]  # Mt is the turnover mass of molecular cores

        beta1 = 1 - alpha1

        if maxV == None:
            # normalFactor1=beta1/( 0 -minV**beta1 )
            normalFactor1 = beta1 / (-minV ** beta1)

        else:
            normalFactor1 = (alpha1 - 1) * minV ** (alpha1 - 1)  # beta1/( maxV**beta1 -minV**beta1 )

        return len(dataArraym) * np.log(normalFactor1) - alpha1 * np.sum(np.log(dataArraym))



    def doIntSingle(self, y, yError, alpha, minV,maxV ):
        """
        minV,and maxV is the integration range

        :param normalFactor:
        :param alpha:
        :param y:
        :param yError:
        :return:
        """

        intMinV = np.max([minV, y - 5 * yError])

        # intMaxV = np.min([maxV,y+5*yError  ])

        intMaxV = y + 5 * yError



        num=np.linspace(intMinV,intMaxV,num=500 )
        #dx=num[1]-num[0]

        functionValues=  num**(-alpha)*np.exp( - (y-num)**2*0.5/yError**2 )

        intAnother= cumtrapz(functionValues,num,initial=1.0)


        if   intAnother[-1] <0:

            print "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW"
        return  intAnother[-1]




    def logSumOfInt(self,yList,yErrorList,alpha,minV,maxV):
        """

        :param yList:
        :param yErrorList:
        :param alpha:
        :param min:
        :param maxV:
        :return:
        """

        sumV=0

        for yS,yErrorS in zip(yList,yErrorList):

            sumV=sumV+np.log( self.doIntSingle(yS,yErrorS,alpha,minV,maxV))


        return sumV


    # fitting function
    def calPowerLawProbLogcomponent1WithError(self ,theta, dataArraym ,dataArrayError , minV,maxV  ):
        """
        This function take the error measure data into account, which is only used to
        #
        :param theta:
        :param dataArraym:
        :param dadtaArrayError:
        :param minV:
        :param maxV:
        :return:
        """

        alpha1 = theta[0]  # Mt is the turnover mass of molecular cores

        beta1 = 1 - alpha1

        if maxV == None:
            # normalFactor1=beta1/( 0 -minV**beta1 )
            normalFactor1 = beta1 / (-minV ** beta1)

        else:
            normalFactor1 = (alpha1 - 1) * minV ** (alpha1 - 1)  # beta1/( maxV**beta1 -minV**beta1 )



        return len(dataArraym) * (  np.log(normalFactor1)  -0.5*np.log(2*np.pi)  )- np.sum(np.log(dataArrayError)) +self.logSumOfInt(dataArraym,dataArrayError,alpha1,minV,maxV )   #- alpha1 * np.sum(np.log(dataArraym))



    def fitPowerLawWithMCMCcomponent1WithError(self, dataArray, dataArrayError, processID,returnDict, sampleN=1000, burn_in=100, minV=None, maxV=None, thin=15):
        """
        fit a power law with MCMC
        """
        # if minV==None or maxV==None:
        # minV= min(dataArray)
        # maxV= max(dataArray)

        mass  = dataArray

        #print "The minimum and maximum masses are (in solar mass), in the MWISPcloud folder", minV, maxV

        np.random.seed()

        alpha1 = np.random.uniform(1.2, 4)  # np.random.exponential(1)

        theta = [alpha1]

        p0 = self.calPowerLawProbLogcomponent1WithError(theta, mass,dataArrayError, minV, maxV)

        sampleK = []
        sampleAlpha = [alpha1]
        widgets = ['MCMCSmapleSlope: ', Percentage(), ' ', Bar(marker='>', left='|', right='|'),
                   ' ', ETA(), ' ', FileTransferSpeed()]  # see docs for other options

        pbar = ProgressBar(widgets=widgets, maxval=sampleN + burn_in + 1)
        pbar.start()

        recordSample = []

        for i in range(100000):

            newAlpha1 = sampleAlpha[-1] + np.random.normal(0, 0.5)  # np.random.exponential(1)

            if newAlpha1<1:
                continue


            theta = [newAlpha1]

            p1 = self.calPowerLawProbLogcomponent1WithError(theta, mass,dataArrayError, minV, maxV)



            if np.isinf(p1):
                continue

            randomR = np.random.uniform(0, 1)

            if p1 >= p0 or p1 - p0 > np.log(randomR):
                p0 = p1

                sampleAlpha.append(newAlpha1)



            else:
                sampleAlpha.append(sampleAlpha[-1])

            if i % thin == 0:


                #print sampleAlpha[-1],"???????????"
                recordSample.append(sampleAlpha[-1])

            pbar.update(len(recordSample))  # this adds a little symbol at each iteration

            if len(recordSample) >= sampleN + burn_in:
                break
        pbar.finish()
        # print mean( sampleAlpha[burn_in:] ), np.std(sampleAlpha[burn_in:]  )

        returnDict[processID] = np.array(recordSample[burn_in:])

        return  #np.array(recordSample[burn_in:])

        # sample theta

    # sample theta

    def fitPowerLawWithMCMCcomponent1(self, dataArray, sampleN=1000, burn_in=100, minV=None, maxV=None, thin=15):
        """
        fit a power law with MCMC
        """
        # if minV==None or maxV==None:
        # minV= min(dataArray)
        # maxV= max(dataArray)

        mass  = dataArray

        print "The minimum and maximum masses are (in solar mass), in the MWISPcloud folder", minV, maxV

        np.random.seed()

        alpha1 = np.random.uniform(1, 5)  # np.random.exponential(1)

        theta = [alpha1]

        p0 = self.calPowerLawProbLogcomponent1(theta, mass, minV, maxV)

        sampleK = []
        sampleAlpha = [alpha1]
        widgets = ['MCMCSmapleSlope: ', Percentage(), ' ', Bar(marker='>', left='|', right='|'),
                   ' ', ETA(), ' ', FileTransferSpeed()]  # see docs for other options

        pbar = ProgressBar(widgets=widgets, maxval=sampleN + burn_in + 1)
        pbar.start()

        recordSample = []

        for i in range(100000):

            newAlpha1 = sampleAlpha[-1] + np.random.normal(0, 0.5)  # np.random.exponential(1)

            theta = [newAlpha1]

            p1 = self.calPowerLawProbLogcomponent1(theta, mass, minV, maxV)

            if np.isinf(p1):
                continue

            randomR = np.random.uniform(0, 1)

            if p1 >= p0 or p1 - p0 > np.log(randomR):
                p0 = p1;

                sampleAlpha.append(newAlpha1)



            else:
                sampleAlpha.append(sampleAlpha[-1])

            if i % thin == 0:
                recordSample.append(sampleAlpha[-1])

            pbar.update(len(recordSample))  # this adds a little symbol at each iteration

            if len(recordSample) > sampleN + burn_in:
                break
        pbar.finish()
        # print mean( sampleAlpha[burn_in:] ), np.std(sampleAlpha[burn_in:]  )
        return np.array(recordSample[burn_in:])

        # sample theta

    # sample theta




    def getAlphaWithMCMCWithErrorMultiChains(self, areaArray, areaArrayData,  minArea=0.03, maxArea=None, sampleN=500, verbose=True, plotTest=False,    saveMark=""):


        print "Fitting index with MCMC and errors are included... {} chains".format(self.nChain )


        if maxArea != None:
            select = np.logical_and(areaArray >= minArea, areaArray <= maxArea)

        else:
            select = areaArray >= minArea

        rawArea = areaArray[select]
        rawAreaError = areaArrayData[select]

        print "Calculating {} chains and each chain has {} samples, thinned by {}.".format(  self.nChain ,   self.chainSample,  self.thinning )

        procs = []

        manager = multiprocessing.Manager()
        returnSampleDic = manager.dict()

        # instantiating process with arguments
        for name in range(self.nChain):
            # print(name)
            #print "Starting process ", name
            # sampleN=1000, burn_in=100, minV=None, maxV=None, thin=15):
            proc = multiprocessing.Process(target=self.fitPowerLawWithMCMCcomponent1WithError, args=(  rawArea, rawAreaError,name, returnSampleDic, self.chainSample,self.burn_in,minArea,maxArea,self.thinning  ) )
            procs.append(proc)
            proc.start()

        for proc in procs:
            proc.join()

        combineSampleArray = []

        for eachArray in returnSampleDic:


            combineSampleArray.append(    returnSampleDic[eachArray] )

        combineSampleArray = np.concatenate(combineSampleArray)

        print len( combineSampleArray),np.mean( combineSampleArray), np.std( combineSampleArray,ddof=1)


        return np.mean( combineSampleArray), np.std( combineSampleArray,ddof=1)

    def getAlphaWithMCMCWithError(self, areaArray, areaArrayData, minArea=0.03, maxArea=None, sampleN=500, verbose=True, plotTest=False,    saveMark=""):
        """
        the area areaArray is acctually the physcial qualities
        :param areaArray:
        :param minArea:
        :param maxArea:
        :return:
        """


        if verbose:
            print "Run first chain for {} molecular clouds.".format(len(rawArea))
        part1 = self.fitPowerLawWithMCMCcomponent1WithError(rawArea, rawAreaError, sampleN=sampleN,  minV=minArea, maxV=maxArea)
        if verbose:
            print "Run second chain for {} molecular clouds.".format(len(rawArea))

        part2 = self.fitPowerLawWithMCMCcomponent1WithError(rawArea, rawAreaError, sampleN=sampleN, minV=minArea, maxV=maxArea)

        allSample = np.concatenate([part1, part2])

        # test plot
        if plotTest:
            fig = plt.figure(figsize=(12, 6))
            ax0 = fig.add_subplot(1, 1, 1)
            # fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
            rc('text', usetex=True)
            rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})

            ax0.scatter(part1, part2, s=10)

            plt.savefig("mcmcSampleTest.pdf", bbox_inches='tight')
            aaaaaa

        meanAlpha = np.mean(allSample)
        stdAlpha = np.std(allSample, ddof=1)
        if verbose:
            print "Alpha Mean: {:.2f}; std: {:.2f}".format(meanAlpha, stdAlpha)

        return round(meanAlpha, 2), round(stdAlpha, 2)


    def getAlphaWithMCMC(self, areaArray, minArea=0.03, maxArea=None, physicalArea=True, verbose=True, plotTest=False,
                         saveMark=""):
        """
        areaArray should be in square armin**2
        :param areaArray:
        :param minArea:
        :param maxArea:
        :return:
        """

        print "Fitting index with MCMC..."

        if not physicalArea:
            areaArray = areaArray / 3600.

        if maxArea != None:
            select = np.logical_and(areaArray >= minArea, areaArray <= maxArea)

        else:
            select = areaArray > minArea

        rawArea = areaArray[select]

        if verbose:
            print "Run first chain for {} molecular clouds.".format(len(rawArea))
        part1 = self.fitPowerLawWithMCMCcomponent1(rawArea, minV=minArea, maxV=maxArea)
        if verbose:
            print "Run second chain for {} molecular clouds.".format(len(rawArea))

        part2 = self.fitPowerLawWithMCMCcomponent1(rawArea, minV=minArea, maxV=maxArea)

        allSample = np.concatenate([part1, part2])

        # test plot
        if plotTest:
            fig = plt.figure(figsize=(12, 6))
            ax0 = fig.add_subplot(1, 1, 1)
            # fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
            rc('text', usetex=True)
            rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})

            ax0.scatter(part1, part2, s=10)

            plt.savefig("mcmcSampleTest.pdf", bbox_inches='tight')
            aaaaaa

        meanAlpha = np.mean(allSample)
        stdAlpha = np.std(allSample, ddof=1)
        if verbose:
            print "Alpha Mean: {:.2f}; std: {:.2f}".format(meanAlpha, stdAlpha)

        return round(meanAlpha, 2), round(stdAlpha, 2)


#test
if 0:
    def rndm(a, b, g, size=1):
        """Power-law gen for pdf(x)\propto x^{g-1} for a<=x<=b"""
        r = np.random.random(size=size)
        ag, bg = a ** g, b ** g
        return (ag + (bg - ag) * r) ** (1. / g)

    #samples=np.random.power(-1.5)

    doPowerLaw=powerLaw1()

    samples= rndm(1,1000,-3.5,2000)

    doPowerLaw.getAlphaWithMCMC(samples,minArea=1,maxArea=100  )
