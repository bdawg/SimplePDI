import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt
from copy import copy
import polzimtools as plz
import time


######################################################################################################
# Settings for reading the raw FITS files (cubed from IDL)
fileExtn = '.fits'

dataPath = '/import/silo4/snert/VAMPIRESData_201607/Analysis/HD169142_21072016/'
filePref = 'cube_HD169142_01_20160721_675-50_EmptySlot_'

#dataPath = '/Volumes/silo4/snert/VAMPIRESData201607/Analysis/HD317501--cal--_21072016/'
#filePref = 'cube_HD317501_01_20160721_675-50_EmptySlot_'

dataPath = '/import/silo4/snert/VAMPIRESData_201609/Analysis/ABAur/'
filePref = 'cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_'
nFiles = 164 *4 # NB actual files, not filenums. So 4x num raw files

dataPath = '/import/silo4/snert/VAMPIRESData_201609/Analysis/MWC758/'
filePref = 'cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_'
nFiles = 172 *4 # NB actual files, not filenums. So 4x num raw files

cubeInfoFile = 'cubeinfoOct2016.idlvar'

saveFilename = 'allSummedImsCubeB_'+filePref

# Assumes HWP angles go as 0, 22.5, 45, 67.5
startFileNum = 0
nFiles = 164 *4 # NB actual files, not filenums. So 4x num raw files
nFiles = 172 *4 # NB actual files, not filenums. So 4x num raw files

showAllImages = False
showPlots = True

centroidThresh = 5000 #Values below this clipped to zero for centroiding
#centroidThresh = 2000

######################################################################################################
# Settings for making the differential images

#bsRot = 74.74 #Differential rotation of beamsplitter, in degrees
bsRot = 89.93 #Why is it not 74.74, like in make_tmpl????????


doMask = False
maskRad = 10 #Radius of central mask (px)
maskVal = 1  #Value to set central masked region to

clim=[-5, 5]
clim=[-100, 100]

dataPath = "/Users/bnorris/DontBackup/simplePDIdata/"

# Make the various difference and ratio images
# loadFilename = 'allSummedImsCube_176files__cube_HD169142_01_20160721_675-50_EmptySlot_.npy'
# loadFilename = 'allSummedImsCube2_cube_HD169142_01_20160721_675-50_EmptySlot_CentroidLLim500.npy'
# loadFilename = 'allSummedImsCube_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
#loadFilename = 'allSummedImsCube_cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_.npy'
#loadFilename = 'allSummedImsCube_cube_HD19820_01_20160918_750-50_EmptySlot_.npz'
loadFilename = 'allSummedImsCubeB_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
loadFilename = 'allSummedImsCube_maxMeth_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
#loadFilename = 'allSummedImsCube_centMeth_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
# loadFilename = 'allSummedImsCube_cenMeth_cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_.npz'
# loadFilename = 'allSummedImsCube_maxMeth_cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_.npz'
#loadFilename = dataPath+loadFilename
#loadFilename = 'allSummedImsCube_maxMeth_cube_AUMic20161015_01_Combined_20161015_750-50_EmptySlot_.npz'

doDerotate = True

# Only use a subset of images? This is in *HWP set* numbering
imRange = [5, 41]

# Reload cube if allSummedIms already exists?
#reload = False


######################################################################################################


class inputImages:
    def __init__(self):
        self.allSummedImsOrig = None
        self.allPAsOrig = None
        self.allStokesIms = None
        self.allPolzstateIms = None


    def loadCube(self, dataPath, loadFilename):
        npzfile = np.load(dataPath+loadFilename)
        self.allSummedImsOrig = npzfile['arr_0']
        self.allPAsOrig = npzfile['arr_1']
        del (npzfile)
        print 'Cube loaded.'


    def processRawFiles(self, dataPath, filePref, cubeInfoFile, nFiles, startFileNum = 0,
                          fileExtn = '.fits', saveFilePref = 'allSummedImsCube_', showImage = False,
                          centroidThresh=5000, saveCube = True, method = 'com'):
        # Read a FITS file to get sizes
        curFilenumStr = '%05d' % startFileNum
        curBSStr = '_1'
        curLCVRStr = '_A'
        curFilename = dataPath + filePref + curFilenumStr + curBSStr + curLCVRStr \
                      + fileExtn
        hdulist = fits.open(curFilename)
        curHDU = hdulist[0]
        curCube = np.transpose(curHDU.data)
        nFrms = curCube.shape[2]
        dim = curCube.shape[0]

        # Get PAs
        cubeInfoFileString = dataPath + cubeInfoFile
        cubeInfo = self.readCubeInfo(cubeInfoFileString)
        pas = cubeInfo.pa

        # Indexes are [:, :, hwpSet, HWP, Channel, LCVRstate]
        # So each set of 4 VAMPIRES on-sky files corresponds to 1 hwpSet.

        allSummedIms = np.zeros([dim, dim, nFiles / 16, 4, 2, 2])
        #curFileNum = startFileNum
        hwpState = 0  # This counts through 4 positions
        curHWPSet = 0  # Increments for each new set of HWP posns

        for f in range(0, nFiles / 4):
            curFileNum = f

            for c in range(0, 2):
                curBSChan = c

                for l in range(0, 2):
                    curLCVRState = l

                    # Generate current filename
                    if curBSChan == 0:
                        curBSStr = '_1'
                    else:
                        curBSStr = '_2'

                    if curLCVRState == 0:
                        curLCVRStr = '_A'
                    else:
                        curLCVRStr = '_B'

                    curFilenumStr = '%05d' % curFileNum
                    curFilename = dataPath + filePref + curFilenumStr + curBSStr + curLCVRStr \
                                  + fileExtn
                    print 'Reading file %s' % curFilename
                    print 'HWPSet %d -- HWP State %d, Channel %d, LCVR %d' % \
                          (curHWPSet, hwpState, curBSChan, curLCVRState)
                    hdulist = fits.open(curFilename)
                    curHDU = hdulist[0]
                    curCube = np.transpose(curHDU.data)

                    # Show summed image
                    curCubeSummed = np.mean(curCube, axis=2)

                    print 'Max val in summed UNshifted image: %f' % np.max(curCubeSummed.ravel())

                    # Make cube of shifted images
                    curCubeShifted = np.zeros([dim, dim, nFrms])
                    # Centre the image
                    for i in range(0, nFrms):

                        llim = centroidThresh
                        curIm = copy(curCube[:, :, i])
                        # print np.min(curIm.min())
                        if method is 'com':
                            curIm[np.where(curIm <= llim)] = 0
                            curCent = ndimage.center_of_mass(curIm)
                        elif method is 'max':
                            curCent = np.where(curIm == curIm.max())
                        else:
                            print 'Error: Unknown centering method specified'
                        # print curCent

                        if showImage:
                            plt.clf()
                            plt.subplot(121)
                            plt.imshow(curCube[:, :, i], interpolation='nearest')
                            plt.plot(curCent[1], curCent[0], 'rx', ms=20, mew=1)
                            plt.subplot(122)
                            plt.imshow(curIm, interpolation='nearest')
                            plt.plot(curCent[1], curCent[0], 'rx', ms=20, mew=1)
                            plt.pause(0.001)

                        imCent = dim / 2
                        rowOffset = imCent - curCent[0]
                        colOffset = imCent - curCent[1]

                        curCubeShifted[:, :, i] = ndimage.interpolation.shift(curCube[:, :, i],
                                                                              [rowOffset, colOffset])

                    curCubeShiftedSummed = np.mean(curCubeShifted, axis=2)
                    curHDU = fits.PrimaryHDU()
                    curHDU.data = curCubeShiftedSummed
                    try:
			curHDU.writeto('autosave.fits', clobber=True)
		    except:
			print "Warning: Couldn't autosave"

                    if showPlots:
                        plt.clf()
                        plt.subplot(121)
                        plt.imshow(curCubeSummed, interpolation='nearest')
                        plt.title('Unshifted sum')
                        plt.subplot(122)
                        plt.title('Shifted sum')
                        plt.imshow(curCubeShiftedSummed, interpolation='nearest')
                        plt.pause(0.001)
                    print 'Max val in summed shifted image:   %f' % np.max(curCubeShiftedSummed.ravel())
                    print ' '
                    allSummedIms[:, :, curHWPSet, hwpState, c, l] = curCubeShiftedSummed

            # Increment HWP state
            hwpState = hwpState + 1
            if hwpState == 4:
                hwpState = 0
                curHWPSet = curHWPSet + 1

        if saveCube:
            # Save the big cube
            saveFilename = saveFilePref + filePref
            np.savez(saveFilename, allSummedIms, pas)

        self.allSummedImsOrig = allSummedIms
        self.allPAsOrig = pas


    class readCubeInfo:
        def __init__(self, cubeInfoFileString):
            # Get useful metadata from cubeinfo file
            cubeinfoObj = io.readsav(cubeInfoFileString, python_dict=False, verbose=False)
            self.UTCs = cubeinfoObj.olog.utc[0]
            self.filters = cubeinfoObj.olog.filter[0]
            self.ras = cubeinfoObj.olog.ra[0]
            self.decs = cubeinfoObj.olog.dec[0]
            self.mask = cubeinfoObj.olog.mask[0]
            self.adate = cubeinfoObj.olog.adate[0]
            self.emgains = cubeinfoObj.olog.emgain[0]
            #self.mffile = cubeinfoObj.plog.mf_file[0]
            self.pkflux = cubeinfoObj.framestats.pkflx[0]
            self.totflux = cubeinfoObj.framestats.totflx[0]
            self.cubename = cubeinfoObj.olog.cube_fname[0][0]
            self.pa = cubeinfoObj.olog.pa[0]
            del (cubeinfoObj)


    def makeDiffIms(self, doMask = False, maskRad = 2., maskVal = 1., clim = [], deRotate = True,
                imRange = [], showPlots = True):
        imSz = self.allSummedImsOrig.shape[0]
        nSets = self.allSummedImsOrig.shape[2]
        if len(imRange) > 0:
            imInd = range(imRange[0], imRange[1])
        else:
            imInd = range(nSets)
        paInds = range(0, len(self.allPAsOrig), 4)
        allPAs = np.asarray(self.allPAsOrig[paInds])
        npa = len(allPAs)
        allPAs = allPAs.reshape((npa/4, 4))
        nUsedSets = len(imInd)

        # Assemble rotated images
        allSummedIms = np.zeros([imSz, imSz, len(imInd), 4, 2, 2])
        pasPerInd = []
        curPaSet = []
        count = 0
        fileCount = 0 #TO keep PAs correct - correspsonds to original file index
        for ind in imInd:
            print "Rotating set: %d" % ind
            for h in range(4):
                pa = allPAs[ind, h]
                curPaSet.append(pa)

                for c in range(2):
                    for l in range(2):
                        curImIn = self.allSummedImsOrig[:, :, ind, h, c, l]
                        # curImIn = allSummedImsOrig[:, :, ind, h, c, l] - h0c0l0med

                        if doMask:
                            maskPosn = imSz / 2
                            yGrid, xGrid = np.ogrid[-maskPosn:imSz - maskPosn, -maskPosn:imSz - maskPosn]
                            maskInds = xGrid ** 2 + yGrid ** 2 <= maskRad ** 2
                            curImIn[maskInds] = maskVal

                        if c == 0:
                            # Rotate channel0 images
                            curImIn = ndimage.interpolation.rotate(np.fliplr(curImIn), bsRot,
                                                                   mode='nearest', reshape=False)
                        if deRotate:
                            allSummedIms[:, :, count, h, c, l] = \
                                ndimage.interpolation.rotate(curImIn, -pa, mode='nearest', reshape=False)
                        else:
                            allSummedIms[:, :, count, h, c, l] = curImIn

            # Approximate set of 4 HWP files as having single PA (the average over files)
            curPa = np.mean(np.asarray(curPaSet))
            curPaSet = []
            pasPerInd.append(curPa)
            count += 1

        pasPerInd = np.asarray(pasPerInd)
        self.pasPerInd = pasPerInd

        # Make a set of Stokes images and derotate their polarisation
        self.allStokesIms = []
        self.allPolzstateIms = []

        for count in range(nUsedSets):
            cur = self.sortPolzState(allSummedIms, count)
            self.allPolzstateIms.append(cur)

            plusQ1 = cur.h0c0l0 - cur.h0c1l0
            minusQ1 = cur.h0c0l1 - cur.h0c1l1
            Q1 = (plusQ1 - minusQ1) / 2
            I0 = (cur.h0c0l0 + cur.h0c1l0 + cur.h0c0l1 + cur.h0c1l1) / 4
            plusQ2 = cur.h2c0l0 - cur.h2c1l0
            minusQ2 = cur.h2c0l1 - cur.h2c1l1
            Q2 = (plusQ2 - minusQ2) / 2
            I2 = (cur.h2c0l0 + cur.h2c1l0 + cur.h2c0l1 + cur.h2c1l1) / 4
            Q = Q1 - Q2
            I_02 = (I0 + I2) / 2
            #pQ = Q / I_02

            plusU1 = cur.h1c0l0 - cur.h1c1l0
            minusU1 = cur.h1c0l1 - cur.h1c1l1
            U1 = (plusU1 - minusU1) / 2
            I1 = (cur.h1c0l0 + cur.h1c1l0 + cur.h1c0l1 + cur.h1c1l1) / 4
            plusU2 = cur.h3c0l0 - cur.h3c1l0
            minusU2 = cur.h3c0l1 - cur.h3c1l1
            U2 = (plusU2 - minusU2) / 2
            I3 = (cur.h3c0l0 + cur.h3c1l0 + cur.h3c0l1 + cur.h3c1l1) / 4
            U = U1 - U2
            I_13 = (I1 + I3) / 2
            #pU = U / I_13
            I = (I_02 + I_13) / 2

            curIm = np.zeros((imSz, imSz, 4))
            curIm[:, :, 0] = I
            curIm[:, :, 1] = Q
            curIm[:, :, 2] = U

            curIm = plz.rotImPolz(curIm, pasPerInd[count])
            self.allStokesIms.append(curIm)

            if showPlots:
                self.plotPolzstates(index=count, fignum=1)
                self.plotStokesIms(index=count, fignum=2)

        allStokesImsArr = np.asarray(self.allStokesIms)
        self.StokesImsSummed = np.mean(allStokesImsArr, axis=0)


    def partialDiffcal(self, calmode):
        """
        Normally triple-differential calibration is performed, via makeDiffIms. This function lets
        you use a subset of polarisation states.
        'calmode' chooses which combination of differential measurements to perform. It uses the
        same numbering system as in diffcal_vampires.
        calType = '0'  ; 0 = Triple Calibration
                ; 1a = Double cal, Woll. + LCVR, HWP0
                ; 1b = Double cal, Woll. + LCVR, HWP45
                ; 2a = Double cal, Woll + HWP, LCVR1
                ; 2b = Double cal, Woll + HWP, LCVR2
                ; 3a = Double cal, LCVR + HWP, Woll.1
                ; 3b = Double cal, LCVR + HWP, Woll.2
                ; 4a = Single cal, Wollaston (LCVR1, HWP0)
                ; 4b = Single cal, Wollaston (LCVR1, HWP45)
                ; 4c = Single cal, Wollaston (LCVR2, HWP0)
                ; 4d = Single cal, Wollaston (LCVR2, HWP45)
                ; 5a = Single cal, LCVR (Woll1, HWP0)
                ; 6a = Single cal, HWP (LCVR1, Woll1)
        """
        allowedCalmodes = ['1a', '2a', '2b', '3a', '4a', '5a', '6a']
        if calmode not in allowedCalmodes:
            print 'Please specify one of the following calmodes:'
            print allowedCalmodes

        imSz = self.allSummedImsOrig.shape[0]
        allPartialIms = []
        count = 0
        for cur in self.allPolzstateIms:
            if calmode is '1a':
                plusQ1 = cur.h0c0l0 - cur.h0c1l0
                minusQ1 = cur.h0c0l1 - cur.h0c1l1
                Q = (plusQ1 - minusQ1) / 2
                plusU1 = cur.h1c0l0 - cur.h1c1l0
                minusU1 = cur.h1c0l1 - cur.h1c1l1
                U = (plusU1 - minusU1) / 2
                Ia = (cur.h0c0l0 + cur.h0c1l0 + cur.h0c0l1 + cur.h0c1l1) / 4
                Ib = (cur.h1c0l0 + cur.h1c1l0 + cur.h1c0l1 + cur.h1c1l1) / 4
                I = (Ia + Ib) / 2
            if calmode is '2a':
                plusQ1 = cur.h0c0l0 - cur.h0c1l0
                minusQ1 = cur.h2c0l0 - cur.h2c1l0
                Q = (plusQ1 - minusQ1) / 2
                plusU1 = cur.h1c0l0 - cur.h1c1l0
                minusU1 = cur.h3c0l0 - cur.h3c1l0
                U = (plusU1 - minusU1) / 2
                Ia = (cur.h0c0l0 + cur.h0c1l0 + cur.h2c0l0 + cur.h2c1l0) / 4
                Ib = (cur.h1c0l0 + cur.h1c1l0 + cur.h3c0l0 + cur.h3c1l0) / 4
                I = (Ia + Ib) / 2
            if calmode is '2b':
                plusQ1 = cur.h0c0l1 - cur.h0c1l1
                minusQ1 = cur.h2c0l1 - cur.h2c1l1
                Q = (plusQ1 - minusQ1) / 2
                plusU1 = cur.h1c0l1 - cur.h1c1l1
                minusU1 = cur.h3c0l1 - cur.h3c1l1
                U = (plusU1 - minusU1) / 2
                Ia = (cur.h0c0l1 + cur.h0c1l1 + cur.h2c0l1 + cur.h2c1l1) / 4
                Ib = (cur.h1c0l1 + cur.h1c1l1 + cur.h3c0l1 + cur.h3c1l1) / 4
                I = (Ia + Ib) / 2
            if calmode is '3a':
                plusQ1 = cur.h0c0l0 - cur.h0c0l1
                minusQ1 = cur.h2c0l0 - cur.h2c0l1
                Q = (plusQ1 - minusQ1) / 2
                plusU1 = cur.h1c0l0 - cur.h1c0l1
                minusU1 = cur.h3c0l0 - cur.h3c0l1
                U = (plusU1 - minusU1) / 2
                Ia = (cur.h0c0l0 + cur.h0c0l1 + cur.h2c0l0 + cur.h2c0l1) / 4
                Ib = (cur.h1c0l0 + cur.h1c0l1 + cur.h3c0l0 + cur.h3c0l1) / 4
                I = (Ia + Ib) / 2
            if calmode is '4a':
                Q = cur.h0c0l0 - cur.h0c1l0
                U = cur.h1c0l0 - cur.h1c1l0
                I = (cur.h0c0l0 + cur.h0c1l0 + cur.h1c0l0 + cur.h1c1l0) / 4
            if calmode is '5a':
                Q = cur.h0c0l0 - cur.h0c0l1
                U = cur.h1c0l0 - cur.h1c0l1
                I = (cur.h0c0l0 + cur.h0c0l1 + cur.h1c0l0 + cur.h1c0l1) / 4
            if calmode is '6a':
                Q = cur.h0c0l0 - cur.h2c0l0
                U = cur.h1c0l0 - cur.h3c0l0
                I = (cur.h0c0l0 + cur.h2c0l0 + cur.h1c0l0 + cur.h3c0l0) / 4

            curIm = np.zeros((imSz, imSz, 4))
            curIm[:, :, 0] = I
            curIm[:, :, 1] = Q
            curIm[:, :, 2] = U
            curIm = plz.rotImPolz(curIm, self.pasPerInd[count])
            allPartialIms.append(curIm)
            count += 1
        allPartialIms = np.asarray(allPartialIms)
        self.PartialImsSummed = np.mean(allPartialIms, axis=0)



    def plotPolzstates(self, data = None, index = None, fignum = 1, clim = None):
        if data == None and index == None:
            print 'Warning: No data or index supplied! Plotting element 0.'
            data = self.allPolzstateIms[0].combined
        elif data == None:
            data = self.allPolzstateIms[index].combined

        idStrings = ['h0c0l0', 'h0c0l1', 'h0c1l0', 'h0c1l1', 'h1c0l0', 'h1c0l1',
                      'h1c1l0', 'h1c1l1', 'h2c0l0', 'h2c0l1', 'h2c1l0', 'h2c1l1',
                      'h3c0l0', 'h3c0l1', 'h3c1l0', 'h3c1l1']

        if showPlots:
            plt.figure(fignum)
            plt.clf()
            for k in range(16):
                im = data[k, :, :]
                plt.subplot(4, 4, k + 1)
                plt.imshow(im, clim=clim, interpolation='Nearest')
                plt.title(idStrings[k])
            plt.pause(0.001)


    def plotStokesIms(self, data = None, index = None, fignum = 1, clim = None):
        if data == None and index == None:
            print 'Warning: No data or index supplied! Plotting element 0.'
            data = self.allStokesIms[0]
        elif data == None:
            data = self.allStokesIms[index]
        pQ = data[:, : ,1] / data[:, :, 0]
        pU = data[:, :, 2] / data[:, :, 0]
        p = np.sqrt(pQ**2 + pU**2)
        idStrings = ['I', 'Q', 'U', 'p', 'pQ', 'pU']

        if showPlots:
            plt.figure(fignum)
            plt.clf()
            for k in range(6):
                if k == 3: im = p
                elif k == 4: im = pQ
                elif k == 5: im = pU
                else:
                    im = data[:, :, k]
                plt.subplot(2, 3, k + 1)
                plt.imshow(im, clim=clim, interpolation='Nearest')
                plt.title(idStrings[k])
            plt.pause(0.001)


    def plotSingleIm(self, data = None, imType = None, fignum = 1, clim = None,
                     radPower = 0, cmap = None, coeffs = None, norm=True, crop=None):
        """
        If no data specified, uses self.StokesImsSummed, and an imType.
        If data and imType are specified, data must be a Stokes image.
        If data and no imType are specified, data must be a flat image.
        """

        if data is None and imType is None:
            print 'Error: You must specify data to plot or an image type.'

        if imType is not None:
            if data is None:
                I = self.StokesImsSummed[:, :, 0]
                Q = self.StokesImsSummed[:, :, 1]
                U = self.StokesImsSummed[:, :, 2]
            else:
                I = data[:, :, 0]
                Q = data[:, :, 1]
                U = data[:, :, 2]

            # By default normalize, but not for fractional pol
            if imType in ['p', 'p2', 'pQ', 'pU']:
                norm = False

            # For toy calibration:
            if coeffs is not None:
                Q = Q * coeffs[0]
                U = U * coeffs[1]

            if imType is 'I':
                data = I
            elif imType is 'Q':
                data = Q
            elif imType is 'U':
                data = U
            elif imType is 'pQ':
                data = Q/I
            elif imType is 'pU':
                data = U/I
            elif imType is 'p':
                data = np.sqrt(Q**2 + U**2) / I
            elif imType is 'p2':
                data = np.sqrt(pQ**2 + pU**2)
            elif imType is 'P':
                data = np.sqrt(Q**2 + U**2)
            else:
                print 'Error: unrecognised image type specified.'

        if radPower != 0:
            sx = data.shape[1]
            sy = data.shape[0]
            inds = np.mgrid[-sy/2:sy/2, -sx/2:sx/2]
            dist = np.sqrt(inds[0,:,:]**2 + inds[1,:,:]**2) + 1
            radMask = dist**radPower
            im = data * radMask
            if norm:
                im /= np.max(im)
        else:
            if norm:
                im = data / np.max(data)
            else:
                im = copy(data)

        if crop is not None:
            sx = data.shape[1]
            sxw = sx*crop / 2
            sy = data.shape[0]
            xi = (sx/2-sxw, sx/2+sxw)
            syw = sy * crop / 2
            yi = (sy / 2 - syw, sy / 2 + syw)
            im = im[yi[0]:yi[1], xi[0]:xi[1]]

        plt.figure(fignum)
        plt.clf()
        plt.imshow(im, clim=clim, interpolation='Nearest', cmap=cmap)
        print 'Image min: %f, max: %f' % (im.min(), im.max())





    class sortPolzState:
        def __init__(self, dataIn, ind):
            self.h0c0l0 = dataIn[:, :, ind, 0, 0, 0]
            self.h0c0l1 = dataIn[:, :, ind, 0, 0, 1]
            self.h0c1l0 = dataIn[:, :, ind, 0, 1, 0]
            self.h0c1l1 = dataIn[:, :, ind, 0, 1, 1]
            self.h1c0l0 = dataIn[:, :, ind, 1, 0, 0]
            self.h1c0l1 = dataIn[:, :, ind, 1, 0, 1]
            self.h1c1l0 = dataIn[:, :, ind, 1, 1, 0]
            self.h1c1l1 = dataIn[:, :, ind, 1, 1, 1]
            self.h2c0l0 = dataIn[:, :, ind, 2, 0, 0]
            self.h2c0l1 = dataIn[:, :, ind, 2, 0, 1]
            self.h2c1l0 = dataIn[:, :, ind, 2, 1, 0]
            self.h2c1l1 = dataIn[:, :, ind, 2, 1, 1]
            self.h3c0l0 = dataIn[:, :, ind, 3, 0, 0]
            self.h3c0l1 = dataIn[:, :, ind, 3, 0, 1]
            self.h3c1l0 = dataIn[:, :, ind, 3, 1, 0]
            self.h3c1l1 = dataIn[:, :, ind, 3, 1, 1]
            self.combined = np.asarray(
                [self.h0c0l0, self.h0c0l1, self.h0c1l0, self.h0c1l1, self.h1c0l0,
                 self.h1c0l1, self.h1c1l0, self.h1c1l1, self.h2c0l0, self.h2c0l1,
                 self.h2c1l0, self.h2c1l1, self.h3c0l0, self.h3c0l1, self.h3c1l0, self.h3c1l1])
            self.idStrings = ['h0c0l0', 'h0c0l1', 'h0c1l0', 'h0c1l1', 'h1c0l0', 'h1c0l1',
                              'h1c1l0', 'h1c1l1', 'h2c0l0', 'h2c0l1', 'h2c1l0', 'h2c1l1',
                              'h3c0l0', 'h3c0l1', 'h3c1l0', 'h3c1l1']






#
p = inputImages()

p.loadCube(dataPath, loadFilename)
p.makeDiffIms(imRange=imRange, showPlots=False)
#p.makeDiffIms(showPlots=False)

I = p.StokesImsSummed[:,:, 0]
Q = p.StokesImsSummed[:,:, 1]
U = p.StokesImsSummed[:,:, 2]
pQ = Q / I
pU = U / I
p2_frac = np.sqrt(pQ**2 + pU**2)
p_frac = np.sqrt(Q**2 + U**2) / I
P = np.sqrt(Q**2 + U**2)

#p.plotSingleIm(imType='p', clim=[-0.05, 0.4])

im = p2_frac
im = np.clip(im,-0.1,0.4)
curHDU = fits.PrimaryHDU()
curHDU.data = im
curHDU.writeto('out.fits', clobber=True)

plt.imshow(im)
plt.pause(0.001)

"""
dataPath = '/import/silo4/snert/VAMPIRESData_201609/Analysis/MWC758/'
filePref = 'cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_'
nFiles = 172 *4 # NB actual files, not filenums. So 4x num raw files

import simplePDIFull_01 as s
cubeInfoFile = 'cubeinfoOct2016.idlvar'
p = s.inputImages()
p.processRawFiles(dataPath, filePref, cubeInfoFile, nFiles,saveFilePref =
    'allSummedImsCube_maxMeth_', method='max')
"""

print 'Done.'







