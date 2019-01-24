"""
VAMPIRES PDI tools
-- This version is to be imported into separate processing scripts, hard-coded tasks removed --
"""


import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt
from copy import copy
import polzimtools as plz
import time





######################################################################################################
# # Usage for reading *raw* files (no IDL pre-cubing)
# ########## NOT YET IMPLEMENTED #############
# dataPath = '/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201709/20170912/'
# fileExtn = '.fits'
# filePref = 'omiCet_01_20170912_750-50_EmptySlot_0'
# startFileBum = 0
# nRawFiles = 4 #68
# nSubFiles = nRawFiles * 2 # To maintain compatibility.
#
# IDLPreCubed = False
# twoCamMode = False
# firstGoodFrame = 2 # Discard frames in a cube before this. For VAMPIRES 2-cam, this should be 2.
#
#
# # General settings:
# showAllImages = False
# showPlots = True
# saveAllSummedIms = False
# centroidThresh = 5000

######################################################################################################


class pdiImages:
    def __init__(self, outDir='./', showPlots=True):
        self.allSummedImsOrig = None
        self.allPAsOrig = None
        self.allStokesIms = None
        self.allPolzstateIms = None
        self.showPlots = showPlots
        self.outDir = outDir
        self.lastRawCube = None


    def loadCube(self, dataPath, loadFilename, twoCamMode = False):
        npzfile = np.load(dataPath+loadFilename)
        self.allSummedImsOrig = npzfile['arr_0']
        self.allPAsOrig = npzfile['arr_1']
        del (npzfile)
        print(self.allPAsOrig)

        if twoCamMode:
            print('loadCube in twoCamMode')
            # TESTING swapping channels for two cam mode

            # # Swap FLC: - maybe this???
            # cur = copy(self.allSummedImsOrig)
            # cur[:, :, :, :, :, 0] = self.allSummedImsOrig[:, :, :, :, :, 1]
            # cur[:, :, :, :, :, 1] = self.allSummedImsOrig[:, :, :, :, :, 0]
            # self.allSummedImsOrig = cur

            # # Swap Channel:
            # cur = copy(self.allSummedImsOrig)
            # cur[:, :, :, :, 0, :] = self.allSummedImsOrig[:, :, :, :, 1, :]
            # cur[:, :, :, :, 1, :] = self.allSummedImsOrig[:, :, :, :, 0, :]
            # self.allSummedImsOrig = cur

            # # Swap HWP:
            # cur = copy(self.allSummedImsOrig)
            # cur[:, :, :, 0, :, :] = self.allSummedImsOrig[:, :, :, 2, :, :]
            # cur[:, :, :, 2, :, :] = self.allSummedImsOrig[:, :, :, 0, :, :]
            # cur[:, :, :, 1, :, :] = self.allSummedImsOrig[:, :, :, 3, :, :]
            # cur[:, :, :, 3, :, :] = self.allSummedImsOrig[:, :, :, 1, :, :]
            # self.allSummedImsOrig = cur

            # # Negate HWP:
            # cur = copy(self.allSummedImsOrig)
            # cur[:, :, :, 1, :, :] = self.allSummedImsOrig[:, :, :, 1, :, :]*-1
            # cur[:, :, :, 3, :, :] = self.allSummedImsOrig[:, :, :, 3, :, :]*-1
            # self.allSummedImsOrig = cur

        print('Cube loaded.')


    def processRawFiles_IDLPrecubed(self, dataPath, filePref, cubeInfoFile, nSubFiles, startFileNum = 0,
                          fileExtn = '.fits', saveFilePref = 'allSummedImsCube_', showImage = False,
                          centroidThresh=5000, saveCube = True, method = 'com',
                          twoCamMode = False, savePath=None):
        # Case when files are sub-cubes made from IDL pipeline
        # TODO - Make show plots behaviout match processRawFiles(), i.e. with self.showPlots
        # TODO - Implement comRad
        if savePath is None:
            savePath = self.outDir

        # Read a FITS file to get sizes
        curFilenumStr = '%05d' % startFileNum
        curBSStr = '_1'
        curLCVRStr = '_A'
        curFilename = dataPath + filePref + curFilenumStr + curBSStr + curLCVRStr \
                      + fileExtn
        hdulist = fits.open(curFilename)
        curHDU = hdulist[0]
        curCube = np.transpose(curHDU.data)
        hdulist.close()
        nFrms = curCube.shape[2]
        dim = curCube.shape[0]

        # Get PAs
        # Note on PAs: As of Oct2017, the PA keyword in vampires data is just the IMR.PAD value,
        # i.e. without the
        # pupil offset (IMR.PAP). To make the PA agree with iObserve values, it *seems* it is:
        # realPA = pad - pap - 180 for Northern targets, or
        # realPA = pad - pap + 180 for Southern targets, or
        cubeInfoFileString = dataPath + cubeInfoFile
        self.cubeInfo = self.readCubeInfo(cubeInfoFileString)
        pas = self.cubeInfo.pa
        # Indexes are [:, :, hwpSet, HWP, Channel, LCVRstate]
        # So each set of 4 VAMPIRES on-sky files corresponds to 1 hwpSet.

        # else:
        #     # Case when files are raw camera files
        #     curFilenumStr = '%d' % startFileNum
        #     curCamStr = '_cam1'
        #     curFilename = dataPath + filePref + curFilenumStr + curCamStr + fileExtn
        #     hdulist = fits.open(curFilename)
        #     curHDU = hdulist[0]
        #     curCube = np.transpose(curHDU.data)
        #     pa = hdulist[1].header['PAD'] #TODO - Do proper PA rather than just PAD value
        #     hdulist.close()
        #     nFrms = curCube.shape[2]
        #     dim = curCube.shape[0]


        allSummedIms = np.zeros([dim, dim, nSubFiles / 16, 4, 2, 2])
        #curFileNum = startFileNum
        hwpState = 0  # This counts through 4 positions
        curHWPSet = 0  # Increments for each new set of HWP posns

        for f in range(0, nSubFiles / 4):
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
                    print('Reading file %s' % curFilename)
                    print('HWPSet %d -- HWP State %d, Channel %d, LCVR %d' % \
                          (curHWPSet, hwpState, curBSChan, curLCVRState))
                    hdulist = fits.open(curFilename)
                    curHDU = hdulist[0]
                    curCube = np.transpose(curHDU.data)
                    hdulist.close()

                    # Show summed image
                    curCubeSummed = np.mean(curCube, axis=2)

                    print('Max val in summed UNshifted image: %f' % np.max(curCubeSummed.ravel()))

                    # Make cube of shifted images
                    curCubeShifted = np.zeros([dim, dim, nFrms])
                    # Centre the image
                    for i in range(0, nFrms):

                        llim = centroidThresh
                        curIm = copy(curCube[:, :, i])
                        # print(np.min(curIm.min()))
                        if method is 'com':
                            curIm[np.where(curIm <= llim)] = 0
                            curCent = ndimage.center_of_mass(curIm)
                        elif method is 'max':
                            curCent = np.where(curIm == curIm.max())
                        else:
                            print('Error: Unknown centering method specified')
                        # print(curCent)

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
                        curHDU.writeto(savePath+'autosave.fits', overwrite=True)
                    except:
                        print("Warning: Couldn't autosave")

                    if self.showPlots:
                        plt.clf()
                        plt.subplot(121)
                        plt.imshow(curCubeSummed, interpolation='nearest')
                        plt.title('Unshifted sum')
                        plt.subplot(122)
                        plt.title('Shifted sum')
                        plt.imshow(curCubeShiftedSummed, interpolation='nearest')
                        plt.pause(0.001)
                    print('Max val in summed shifted image:   %f' % np.max(curCubeShiftedSummed.ravel()))
                    print(' ')
                    allSummedIms[:, :, curHWPSet, hwpState, c, l] = curCubeShiftedSummed

            # Increment HWP state
            hwpState = hwpState + 1
            if hwpState == 4:
                hwpState = 0
                curHWPSet = curHWPSet + 1

        if saveCube:
            # Save the big cube
            saveFilename = savePath + saveFilePref + filePref
            print('Saving processed cube as '+saveFilename)
            np.savez(saveFilename, allSummedIms, pas)

        self.allSummedImsOrig = allSummedIms
        self.allPAsOrig = pas



    def processRawFiles(self, dataPath, filePref, nSubFiles, startFileNum = 0,
                        fileExtn = '.fits', saveFilePref = 'allSummedImsCube_', showPlots=None,
                        centroidThresh=5000, saveCube = True, method = 'com', comRad=0,
                        twoCamMode = True, savePath=None, bgMode=0, bgVals=(0, 0),
                        darkFilename=None, skyAnn=(0, 0), showAllIms=False, previewCropRad=None,
                        figNums=(1,2,3), luckyCriteria = None,  # 'pkflux', 'totflux' or 'l2flux' or None
                        luckyPercent = 99, singleFileOnly=False):
        # Case when files are original FITS files pipeline
        # Note - currently only works with two-camera data

        # if comRad > 0 then com will only look for COM in a box of this half-width centred at max
        # If previewCropRad is not None then crop preview images to this half-width

        if not twoCamMode:
            print('ERROR: Cannot do single-cam data without IDL pre-cubing.')
            return

        if savePath is None:
            savePath = self.outDir

        if bgMode == 2:
            npzfile = np.load(darkFilename)
            darkFrames = npzfile['finalDarks']
            del (npzfile)

        if showPlots is None:
            showPlots = self.showPlots

        # Read a FITS file to get sizes
        curFilenumStr = '%d' % startFileNum
        curCamStr = '_cam1'
        curFilename = dataPath + filePref + curFilenumStr + curCamStr + fileExtn
        hdulist = fits.open(curFilename)
        curHDU = hdulist[0]
        curCube = np.transpose(curHDU.data)
        nFrms = curCube.shape[2]
        dim = curCube.shape[0]

        # Indexes are [:, :, hwpSet, HWP, Channel (camera), FLCstate]
        # So each set of 4 VAMPIRES on-sky files corresponds to 1 hwpSet.
        if singleFileOnly:
            nSubFiles=16
        allSummedIms = np.zeros([dim, dim, nSubFiles // 16, 4, 2, 2])
        allRawPAs = []  # List of PAs from each file read
        # curFileNum = startFileNum
        hwpState = 0  # This counts through 4 positions
        curHWPSet = 0  # Increments for each new set of HWP posns

        for f in range(0, nSubFiles // 4):
            curFileNum = f

            if singleFileOnly:
                print('WARNING: Only using first file, and duplicating it.')
                curFileNum = 0

            for c in range(0, 2):
                curChan = c

                # Generate current filename
                if curChan == 0:
                    curCamStr = '_cam1'
                else:
                    curCamStr = '_cam2'

                curFilenumStr = '%d' % curFileNum
                curFilename = dataPath + filePref + curFilenumStr + curCamStr + fileExtn

                print('Reading file %s' % curFilename)
                # print 'HWPSet %d -- HWP State %d, Channel %d, LCVR %d' % \
                #       (curHWPSet, hwpState, curChan, curFLCState)
                hdulist = fits.open(curFilename)
                curHDU = hdulist[0]
                curSuperCube = np.transpose(curHDU.data)  # 'Super' since both FLC states

                if twoCamMode & curChan == 1:
                    curSuperCube = np.fliplr(curSuperCube)

                # Get PA from header (use PAD, and assume instrument offset added later)
                try:
                    pa = hdulist[1].header['PAD']
                    allRawPAs.append(pa)
                except:
                    print('WARNING: Could not get PA from header')
                    allRawPAs.append(0)

                # Now sort frames into the 2 FLC states and discard 1st 2 frames
                nSubFrms = nFrms // 2 - 1
                curCube_FLC1 = np.zeros((dim, dim, nSubFrms))
                curCube_FLC2 = np.zeros((dim, dim, nSubFrms))
                curFLC = 1
                count = 0
                for k in range(2, nFrms):
                    if curFLC == 1:
                        curCube_FLC1[:, :, count] = curSuperCube[:, :, k]

                    else:
                        curCube_FLC2[:, :, count] = curSuperCube[:, :, k]
                        count = count + 1

                    if curFLC == 1:
                        curFLC = 2
                    else:
                        curFLC = 1

                for l in range(0, 2):
                    curFLCState = l
                    print('HWPSet %d -- HWP State %d, Channel %d, LCVR %d' % \
                          (curHWPSet, hwpState, curChan, curFLCState))

                    if l == 0:
                        curCube = curCube_FLC1
                    else:
                        curCube = curCube_FLC2

                    # Lucky imaging - pick the best ones
                    # Evaluate based on cam1, and use the same selection for both cams
                    if luckyCriteria is not None:
                        unluckyCubeSummed = np.mean(curCube, axis=2)
                        print('Max val in summed Unlucky image: %f' % np.max(unluckyCubeSummed.ravel()))

                        if curChan == 0:
                            if luckyCriteria is 'pkflux':
                                luckMet = np.amax(curCube,axis=(0,1))
                            elif luckyCriteria is 'totflux':
                                luckMet = np.sum(curCube, axis=(0, 1))
                            elif luckyCriteria is 'l2flux':
                                luckMet = np.sum(curCube**2,axis=(0,1))
                            elif luckyCriteria is None:
                                pass
                            else:
                                print('Error - invalid lucky criteria specified')
                                return

                            cut = np.percentile(luckMet, luckyPercent)
                            goodInds = np.where(luckMet >= cut)

                            if l == 0:
                                goodInds_FLC1 = goodInds
                            else:
                                goodInds_FLC2 = goodInds

                            print('Lucky imaging - using best %4.1f%%, metric has mean %4.2g, sigma %4.2g' %
                                  (luckyPercent, np.mean(luckMet), np.std(luckMet)))

                        if l == 0:
                            curCube = curCube[:, :, goodInds_FLC1[0]]
                            print('Sum of goodInds used: %d' % np.sum(goodInds_FLC1[0]))
                        else:
                            curCube = curCube[:, :, goodInds_FLC2[0]]
                            print('Sum of goodInds used: %d' % np.sum(goodInds_FLC2[0]))
                    else:
                        unluckyCubeSummed = np.zeros((dim, dim))

                    self.lastRawCube = curCube # For debugging

                    curCubeSummed = np.mean(curCube, axis=2)
                    print('Max val in summed UNshifted image: %f' % np.max(curCubeSummed.ravel()))

                    if method is 'none':
                        curCubeShifted = curCube
                    else:
                        # Make cube of shifted images
                        curCubeShifted = np.zeros_like(curCube)
                        # Centre the image
                        for i in range(0, curCube.shape[2]):
                            llim = centroidThresh
                            curIm = copy(curCube[:, :, i])
                            # print np.min(curIm.min())
                            if method is 'com':
                                if comRad == 0:
                                    curIm[np.where(curIm <= llim)] = 0
                                    curCent = ndimage.center_of_mass(curIm)
                                else:
                                    curIm[np.where(curIm <= llim)] = 0
                                    roughCent = np.where(curIm == curIm.max())
                                    mask = np.zeros_like(curIm)
                                    mask[roughCent[0][0] - comRad:roughCent[0][0] + comRad,
                                        roughCent[1][0] - comRad:roughCent[1][0] + comRad] = 1
                                    curImMasked = curIm * mask
                                    curCent = ndimage.center_of_mass(curImMasked)
                            elif method is 'max':
                                curCent = np.where(curIm == curIm.max())
                            elif method is 'none':
                                pass
                            else:
                                print('Error: Unknown centering method specified')
                            # print curCent

                            if showAllIms:
                                plt.clf()
                                plt.subplot(121)
                                im = curCube[:, :, i]
                                plt.imshow(im, interpolation='nearest')
                                plt.plot(curCent[1], curCent[0], 'rx', ms=50, mew=1)
                                if previewCropRad is not None:
                                    cent = curCube.shape[0] // 2
                                    pc = previewCropRad
                                    plt.xlim(cent - pc, cent + pc)
                                    plt.ylim(cent - pc, cent + pc)
                                plt.subplot(122)
                                im = curIm
                                plt.imshow(im, interpolation='nearest')
                                plt.plot(curCent[1], curCent[0], 'rx', ms=50, mew=1)
                                if previewCropRad is not None:
                                    plt.xlim(cent - pc, cent + pc)
                                    plt.ylim(cent - pc, cent + pc)
                                plt.pause(0.001)

                            if bgMode == 2:  # Darkframe subtraction must happen before shifting
                                curCube[:, :, i] = curCube[:, :, i] - darkFrames[:, :, curChan]
                                darkMed = np.median(darkFrames[:, :, curChan])

                            imCent = dim / 2
                            rowOffset = imCent - curCent[0]
                            colOffset = imCent - curCent[1]

                            curCubeShifted[:, :, i] = ndimage.interpolation.shift(curCube[:, :, i],
                                                                                  [rowOffset, colOffset])
                    curCubeShiftedSummed = np.mean(curCubeShifted, axis=2)
                    curHDU = fits.PrimaryHDU()
                    curHDU.data = curCubeShiftedSummed
                    try:
                        curHDU.writeto('autosave.fits', overwrite=True)
                    except:
                        print("Warning: Couldn't autosave")

                    if showPlots:
                        plt.clf()
                        plt.subplot(131)
                        plt.imshow(unluckyCubeSummed, interpolation='nearest')
                        if previewCropRad is not None:
                            cent = curCube.shape[0] // 2
                            pc = previewCropRad
                            plt.xlim(cent - pc, cent + pc)
                            plt.ylim(cent - pc, cent + pc)
                        plt.title('Unlucky, unshifted sum')
                        plt.subplot(132)
                        plt.imshow(curCubeSummed, interpolation='nearest')
                        if previewCropRad is not None:
                            plt.xlim(cent - pc, cent + pc)
                            plt.ylim(cent - pc, cent + pc)
                        plt.title('Unshifted sum')
                        plt.subplot(133)
                        plt.title('Shifted sum')
                        plt.imshow(curCubeShiftedSummed, interpolation='nearest')
                        if previewCropRad is not None:
                            plt.xlim(cent - pc, cent + pc)
                            plt.ylim(cent - pc, cent + pc)
                        plt.pause(0.001)
                    print('Max val in summed shifted image:   %f' % np.max(curCubeShiftedSummed.ravel()))
                    print(' ')

                    if bgMode == 1:
                        curCubeShiftedSummed = curCubeShiftedSummed - bgVals[curChan]
                        print('Dark value subtracted: %f' % bgVals[curChan])

                    if bgMode == 2:
                        print('Dark frame subtracted, with median %f' % darkMed)

                    if bgMode == 3:
                        skyAnnInner = skyAnn[0]
                        skyAnnOuter =  skyAnn[1]
                        px = np.arange(dim) - dim / 2
                        xv, yv = np.meshgrid(px, px)
                        rpix = np.sqrt(xv ** 2 + yv ** 2)
                        skyInds = np.where(np.logical_and(rpix >= skyAnnInner, rpix <= skyAnnOuter))
                        skyVals = curCubeShiftedSummed[skyInds]
                        skyVal = np.median(skyVals)
                        print('Sky value subtracted: %f' % skyVal)
                        curCubeShiftedSummed = curCubeShiftedSummed - skyVal

                    allSummedIms[:, :, curHWPSet, hwpState, c, l] = curCubeShiftedSummed

            # Increment HWP state
            hwpState = hwpState + 1
            if hwpState == 4:
                hwpState = 0
                curHWPSet = curHWPSet + 1

        # This is to maintain compatability with the single-camera data
        pas = np.repeat(allRawPAs, 2)

        if saveCube:
            # Save the big cube
            saveFilename = savePath + saveFilePref + filePref
            print('Saving processed cube as '+saveFilename)
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
            imRange = [], showPlots = None, maskfile = [], calMat_PreRot = [],
            calMat_PostRot = [], twoCamMode = True, doFlats = False, Ibias = 0, darkOffset=None,
            saveAllSummedIms=False, bsRot=89.93, suppressWarnings=True):
        if showPlots is None:
            showPlots = self.showPlots
        if Ibias > 0:
            print('WARNING: Stokes I bias applied')
        if suppressWarnings:
            np.seterr(invalid='ignore')
            np.seterr(divide='ignore')

        imSz = self.allSummedImsOrig.shape[0]
        nSets = self.allSummedImsOrig.shape[2]
        if len(imRange) > 0:
            imInd = range(imRange[0], imRange[1])
        else:
            imInd = range(nSets)
        paInds = range(0, len(self.allPAsOrig), 4)
        allPAs = np.asarray(self.allPAsOrig[paInds])
        npa = len(allPAs)
        allPAs = allPAs.reshape((int(npa/4), 4))
        self.allPAsReshaped = allPAs
        nUsedSets = len(imInd)

        if len(maskfile) > 0:
            maskfileobj = np.load(maskfile)
            mask1 = maskfileobj['mask1']
            mask2 = maskfileobj['mask2']
        else:
            mask1 = np.ones((imSz, imSz))
            mask2 = np.ones((imSz, imSz))

        # Assemble rotated images
        allSummedIms = np.zeros([imSz, imSz, len(imInd), 4, 2, 2])
        pasPerInd = []
        curPaSet = []
        count = 0
        fileCount = 0 #TO keep PAs correct - corresponds to original file index
        for ind in imInd:
            print("Rotating set: %d" % ind)
            for h in range(4):
                pa = allPAs[ind, h]
                curPaSet.append(pa)

                for c in range(2):
                    for l in range(2):
                        curImIn = self.allSummedImsOrig[:, :, ind, h, c, l]
                        # curImIn = allSummedImsOrig[:, :, ind, h, c, l] - h0c0l0med

                        if darkOffset is not None:
                            curImIn = curImIn - darkOffset

                        if c == 0:
                            curImIn = curImIn * mask1
                        if c == 1:
                            curImIn = curImIn * mask2

                        if doMask:
                            maskPosn = imSz / 2
                            yGrid, xGrid = np.ogrid[-maskPosn:imSz - maskPosn, -maskPosn:imSz - maskPosn]
                            maskInds = xGrid ** 2 + yGrid ** 2 <= maskRad ** 2
                            curImIn[maskInds] = maskVal

                        if doFlats:
                            # TODO - hacked in filenames etc right now - fix this.
                            flatNpzFile = np.load('flats750Aug2015.npz')
                            flatC1 = flatNpzFile['flatC1']
                            flatC2 = flatNpzFile['flatC2']
                            if c == 0:
                                curImIn = curImIn / flatC1
                            if c == 1:
                                curImIn = curImIn / flatC2

                        if twoCamMode:
                            if c == 1:
                                pass
                                #curImIn = np.fliplr(curImIn)
                                # Instead, do this beforehand during processRawFiles, for alignment.
                        else:
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

        if saveAllSummedIms:
            fname = 'allSummedImsSaved.npz'
            np.savez(fname, allSummedIms=allSummedIms, pasPerInd=pasPerInd)

        # Make a set of Stokes images and derotate their polarisation
        self.allStokesIms = []
        self.allPolzstateIms = []
        self.allPolzratioQ = []
        self.allPolzratioU = []

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
            ############Q = Q1 + Q2
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
            ##############U = U1 + U2
            I_13 = (I1 + I3) / 2
            #pU = U / I_13
            I = (I_02 + I_13) / 2

            curIm = np.zeros((imSz, imSz, 4))
            curIm[:, :, 0] = I + Ibias # Ibias is a kludge to increase visibility of some features
                # in fractional polz map
            curIm[:, :, 1] = Q
            curIm[:, :, 2] = U

            if len(calMat_PreRot) > 0:
                curIm = plz.matIm(curIm, calMat_PreRot)

            curIm = plz.rotImPolz(curIm, pasPerInd[count])

            if len(calMat_PostRot) > 0:
                curIm = plz.matIm(curIm, calMat_PostRot)

            self.allStokesIms.append(curIm)


            # Now calculate fractional polarisations using the triple-ratio method
            RQ1 = np.sqrt((cur.h0c0l0 / cur.h0c1l0) / (cur.h0c0l1 / cur.h0c1l1))
            RQ2 = np.sqrt((cur.h2c0l0 / cur.h2c1l0) / (cur.h2c0l1 / cur.h2c1l1))
            RQ = np.sqrt(RQ1 / RQ2)
            pQ = (RQ - 1) / (RQ + 1)

            RU1 = np.sqrt((cur.h1c0l0 / cur.h1c1l0) / (cur.h1c0l1 / cur.h1c1l1))
            RU2 = np.sqrt((cur.h3c0l0 / cur.h3c1l0) / (cur.h3c0l1 / cur.h3c1l1))
            RU = np.sqrt(RU1 / RU2)
            pU = (RU - 1) / (RU + 1)
            # cur_p_frac = np.sqrt(pQ**2 + pU**2)

            # Make dummy Stokes vector so can reuse the optimised rotImPolz
            curIm = np.zeros((imSz, imSz, 4))
            curIm[:, :, 1] = pQ
            curIm[:, :, 2] = pU

            curIm = plz.rotImPolz(curIm, pasPerInd[count])
            self.allPolzratioQ.append(curIm[:, :, 1])
            self.allPolzratioU.append(curIm[:, :, 2])

            if showPlots:
                self.plotPolzstates(index=count, fignum=1)
                self.plotStokesIms(index=count, fignum=2)

        allStokesImsArr = np.asarray(self.allStokesIms)
        self.StokesImsSummed = np.mean(allStokesImsArr, axis=0)

        summedpQ = np.asarray(self.allPolzratioQ).mean(axis=0)
        summedpU = np.asarray(self.allPolzratioU).mean(axis=0)
        self.polzRatioImsSummed = np.sqrt(summedpQ**2 + summedpU**2)



    def partialDiffcal(self, calmode, calMat_PreRot = [], showImageType=None,
                       showImageCrop=None, clim=None, apertureRad=None, apertureMean=False):
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
                ; 5b = Single cal, LCVR (Woll1, HWP45)
                ; 6a = Single cal, HWP (LCVR1, Woll1)
        """
        allowedCalmodes = ['1a', '2a', '2b', '3a', '3b', '4a', '5a', '5b', '6a']
        if calmode not in allowedCalmodes:
            print('Please specify one of the following calmodes:')
            print(allowedCalmodes)

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
            if calmode is '1b':
                plusQ1 = cur.h2c0l0 - cur.h2c1l0
                minusQ1 = cur.h2c0l1 - cur.h2c1l1
                Q = (plusQ1 - minusQ1) / 2
                plusU1 = cur.h3c0l0 - cur.h3c1l0
                minusU1 = cur.h3c0l1 - cur.h3c1l1
                U = (plusU1 - minusU1) / 2
                Ia = (cur.h2c0l0 + cur.h2c1l0 + cur.h2c0l1 + cur.h2c1l1) / 4
                Ib = (cur.h3c0l0 + cur.h3c1l0 + cur.h3c0l1 + cur.h3c1l1) / 4
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
            if calmode is '3b':
                plusQ1 = cur.h0c1l0 - cur.h0c1l1
                minusQ1 = cur.h2c1l0 - cur.h2c1l1
                Q = (plusQ1 - minusQ1) / 2
                plusU1 = cur.h1c1l0 - cur.h1c1l1
                minusU1 = cur.h3c1l0 - cur.h3c1l1
                U = (plusU1 - minusU1) / 2
                Ia = (cur.h0c1l0 + cur.h0c1l1 + cur.h2c1l0 + cur.h2c1l1) / 4
                Ib = (cur.h1c1l0 + cur.h1c1l1 + cur.h3c1l0 + cur.h3c1l1) / 4
                I = (Ia + Ib) / 2
            if calmode is '4a':
                Q = cur.h0c0l0 - cur.h0c1l0
                U = cur.h1c0l0 - cur.h1c1l0
                I = (cur.h0c0l0 + cur.h0c1l0 + cur.h1c0l0 + cur.h1c1l0) / 4
            if calmode is '5a':
                Q = cur.h0c0l0 - cur.h0c0l1
                U = cur.h1c0l0 - cur.h1c0l1
                I = (cur.h0c0l0 + cur.h0c0l1 + cur.h1c0l0 + cur.h1c0l1) / 4
            if calmode is '5b':
                Q = cur.h2c0l0 - cur.h2c0l1
                U = cur.h3c0l0 - cur.h3c0l1
                I = (cur.h2c0l0 + cur.h2c0l1 + cur.h3c0l0 + cur.h3c0l1) / 4
            if calmode is '6a':
                Q = cur.h0c0l0 - cur.h2c0l0
                U = cur.h1c0l0 - cur.h3c0l0
                I = (cur.h0c0l0 + cur.h2c0l0 + cur.h1c0l0 + cur.h3c0l0) / 4

            curIm = np.zeros((imSz, imSz, 4))
            curIm[:, :, 0] = I
            curIm[:, :, 1] = Q
            curIm[:, :, 2] = U
            if len(calMat_PreRot) > 0:
                curIm = plz.matIm(curIm, calMat_PreRot)

            curIm = plz.rotImPolz(curIm, self.pasPerInd[count])
            allPartialIms.append(curIm)
            count += 1
        allPartialIms = np.asarray(allPartialIms)
        self.PartialImsSummed = np.mean(allPartialIms, axis=0)
        print('Done.')

        if showImageType is not None:
            self.plotSingleIm(data='partial', imType=showImageType, crop=showImageCrop,
                              clim=clim, apertureRad=apertureRad, apertureMean=apertureMean)



    def plotPolzstates(self, data = None, index = None, fignum = 1, clim = None):
        if data == None and index == None:
            print('Warning: No data or index supplied! Plotting element 0.')
            data = self.allPolzstateIms[0].combined
        elif data == None:
            data = self.allPolzstateIms[index].combined

        idStrings = ['h0c0l0', 'h0c0l1', 'h0c1l0', 'h0c1l1', 'h1c0l0', 'h1c0l1',
                      'h1c1l0', 'h1c1l1', 'h2c0l0', 'h2c0l1', 'h2c1l0', 'h2c1l1',
                      'h3c0l0', 'h3c0l1', 'h3c1l0', 'h3c1l1']

        if self.showPlots:
            plt.figure(fignum)
            plt.clf()
            for k in range(16):
                im = data[k, :, :]
                plt.subplot(4, 4, k + 1)
                plt.imshow(im, clim=clim, interpolation='Nearest')
                plt.title(idStrings[k])
            plt.pause(0.001)


    def plotStokesIms(self, data = None, index = None, fignum = 2, clim = None,
                      crop=None, saveFITS=False, savePath=None, showColorbar=True,
                      savePlotPNG=True, apertureRad=None, returnVec=False,
                      outFilename='outStokesFig.png'):
        if data is None and index is None:
            # print('Warning: No data or index supplied! Plotting element 0.')
            # data = self.allStokesIms[0]
            print('Warning: No data or index supplied! Plotting SUMMED Stokes ims.')
            data = self.StokesImsSummed
        elif data is None:
            data = self.allStokesIms[index]
        elif data is 'partial':
            data = self.PartialImsSummed

        if savePath is None:
            savePath = self.outDir

        pQ = data[:, : ,1] / data[:, :, 0]
        pU = data[:, :, 2] / data[:, :, 0]
        p = np.sqrt(pQ**2 + pU**2)
        idStrings = ['I', 'Q', 'U', 'p', 'pQ', 'pU']

        if saveFITS:
            fitsData = []

        if self.showPlots:
            plt.figure(fignum, figsize=[15,8])
            plt.clf()
            fluxes = []
            for k in range(6):
                if k == 3: im = p
                elif k == 4: im = pQ
                elif k == 5: im = pU
                else:
                    im = data[:, :, k]

                if apertureRad is not None:
                    if k in (3, 4, 5):
                        apertureMean = True
                    else:
                        apertureMean = False
                    flux = self.measureAperture(im, apertureRad, mean=apertureMean)
                    fluxes.append(flux)
                    print('Aperture flux for %s: %f' % (idStrings[k], flux))

                if crop is not None:
                    sx = data.shape[1]
                    sxw = int(sx * crop / 2)
                    sy = data.shape[0]
                    xi = (sx // 2 - sxw, sx // 2 + sxw)
                    syw = int(sy * crop / 2)
                    yi = (sy // 2 - syw, sy // 2 + syw)
                    im = im[yi[0]:yi[1], xi[0]:xi[1]]

                if saveFITS:
                    fitsData.append(im)

                plt.subplot(2, 3, k + 1)
                plt.imshow(im, clim=clim, interpolation='Nearest')
                plt.title(idStrings[k])
                if showColorbar:
                    plt.colorbar()

            plt.tight_layout()
            plt.pause(0.001)

        if apertureRad is not None:
            fluxes = np.asarray(fluxes)
            pQ_s = fluxes[1] / fluxes[0]
            pU_s = fluxes[2] / fluxes[0]
            p_s = np.sqrt(pQ_s ** 2 + pU_s ** 2)
            print('Scalar p: %f' % p_s)
            print('Scalar pQ: %f' % pQ_s)
            print('Scalar pU: %f' % pU_s)
            print(' ')
            # Tab-delimited: I Q U p pQ pU p_s pQ_s pU_s
            polzFluxVec = [fluxes[0], fluxes[1], fluxes[2], fluxes[3], fluxes[4],
                           fluxes[5], p_s, pQ_s, pU_s]

        if savePlotPNG:
            outPath = savePath + outFilename
            plt.savefig(outPath, dpi=100)

        if saveFITS:
            fitsData = np.asarray(fitsData)
            outPath = savePath+'outStokesIms.fits'
            print('Writing FITS file to '+outPath)
            fits.writeto(outPath, fitsData, overwrite=True)

        if returnVec:
            return polzFluxVec



    def plotSingleIm(self, data = None, imType = None, fignum = 1, clim = None,
                     radPower = 0, cmap = None, coeffs = None, norm=True, crop=None,
                     showColorbar=True, apertureRad=None, apertureMean=False, pixScale=6.5):
        """
        If no data specified, uses self.StokesImsSummed, and an imType.
        If data and imType are specified, data must be a Stokes image.
        If data and no imType are specified, data must be a flat image.

        apertureRad: Radius of a centred aperture to sum flux in
        """

        if data is None and imType is None:
            print('Error: You must specify data to plot or an image type.')

        if data is 'partial':
            data = self.PartialImsSummed

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
                pQ = Q / I
                pU = U / I
                data = np.sqrt(pQ**2 + pU**2)
            elif imType is 'P':
                data = np.sqrt(Q**2 + U**2)
            else:
                print('Error: unrecognised image type specified.')

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
            sxw = int(sx*crop / 2)
            sy = data.shape[0]
            xi = (sx//2-sxw, sx//2+sxw)
            syw = int(sy * crop / 2)
            yi = (sy // 2 - syw, sy // 2 + syw)
            im = im[yi[0]:yi[1], xi[0]:xi[1]]

        plt.figure(fignum)
        plt.clf()

        bLen = 100./pixScale
        plt.plot((10, 10+bLen), (10, 10), 'w')
        plt.text(10, 17, '~100 mas', color='w')
        # bLen = 100./pixScale
        # plt.plot((5, 5+bLen), (5, 5), 'w')
        # plt.text(5, 7, '~100 mas', color='w')

        plt.imshow(im, clim=clim, interpolation='Nearest', cmap=cmap)
        if showColorbar:
            plt.colorbar()

        if apertureRad is not None:
            flux = self.measureAperture(im, apertureRad, mean=apertureMean)
            print('Aperture flux: %f' % flux)

        print('Image min: %f, max: %f' % (im.min(), im.max()))
        return im



    def measureAperture(self, im, radius, centre=None, mean=False, showIm=False, showImSize=None):
        """
        Return aperture flux
        :param im: Input image
        :param radius: Radius of measurement aperture
        :param centre: Centre of aperture. If None then use image centre.
        :param mean: Set to True to return teh mean, otherwise return the sum.
        :param showIm: If True, show the masked image
        :return: sum or mean of aperture flux
        """
        h = im.shape[0]
        w = im.shape[1]
        if centre is None:
            centre=[h//2, w//2]
        Y, X = np.ogrid[:h, :w]
        dist = np.sqrt((X-centre[0])**2 + (Y-centre[1])**2)
        mask = dist >= radius
        imMasked = np.ma.masked_array(im, mask=mask)
        if mean:
            result = np.mean(imMasked)
        else:
            result = np.sum(imMasked)

        if showIm:
            if showImSize is not None:
                imToShow = imMasked[centre[0]-showImSize//2:centre[0]+showImSize//2-1,
                           centre[1] - showImSize // 2:centre[1] + showImSize // 2 - 1]
            else:
                imToShow = imMasked

            plt.clf()
            plt.imshow(imToShow)

        return result



    def encircEnergy(self, im, radii=(1, 32), centre=None, mean=False, showIm=False, showImSize=None):
        rads = np.arange(radii[0], radii[1])
        fluxes = []
        for r in rads:
            f = self.measureAperture(im, r, centre=centre, mean=mean, showIm=showIm, showImSize=showImSize)
            fluxes.append(f)
            if showIm:
                plt.pause(0.1)
        fluxes = np.asarray(fluxes)
        plt.figure()
        plt.plot(rads, fluxes)
        return fluxes


    def doMultiCals(self, crop=0.25, apertureRad=20):
        # Convenience function to measure and report multiple calibration types
        allPolzVecs = []
        print('Calmode 0:')
        v = self.plotStokesIms(crop=crop, apertureRad=apertureRad, returnVec=True,
                               fignum=1, outFilename='stokesFig_Cal0')
        allPolzVecs.append(v)

        self.partialDiffcal('3a')
        print('Calmode 3a:')
        v = self.plotStokesIms(data='partial', crop=crop, apertureRad=apertureRad,
                               returnVec=True, fignum=2, outFilename='stokesFig_Cal3a')
        allPolzVecs.append(v)

        self.partialDiffcal('3b')
        print('Calmode 3b:')
        v = self.plotStokesIms(data='partial', crop=crop, apertureRad=apertureRad,
                               returnVec=True, fignum=3, outFilename='stokesFig_Cal3b')
        allPolzVecs.append(v)

        print('Combined vector (0, 3a, 3b):')
        allPolzVecsString = ''
        for k in allPolzVecs:
            for l in k:
                allPolzVecsString = allPolzVecsString + '%f'%l + '\t'
        print(allPolzVecsString)



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
