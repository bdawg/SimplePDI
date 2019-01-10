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

IDLPreCubed = True # True for IDL pipeline cubes, False for raw camera files (only for two-cam data)

# dataPath = '/import/silo4/snert/VAMPIRESData_201607/Analysis/HD169142_21072016/'
# filePref = 'cube_HD169142_01_20160721_675-50_EmptySlot_'

#dataPath = '/Volumes/silo4/snert/VAMPIRESData201607/Analysis/HD317501--cal--_21072016/'
#filePref = 'cube_HD317501_01_20160721_675-50_EmptySlot_'

#dataPath = '/import/silo4/snert/VAMPIRESData_201609/Analysis/ABAur/'
dataPath = "/Users/bnorris/DontBackup/simplePDIdata/ABAur/"
filePref = 'cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_'
nSubFiles = 164 *4 # NB actual files, not filenums. So 4x num raw files
cubeInfoFile = 'cubeinfoOct2016.idlvar'

# dataPath = '/import/silo4/snert/VAMPIRESData_201609/Analysis/MWC758/'
# filePref = 'cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_'
# nSubFiles = 172 *4 # NB actual files, not filenums. So 4x num raw files

# dataPath = "/Users/bnorris/DontBackup/simplePDIdata/ABAur_201610/"
# filePref = 'cube_ABAur_01-02Combined_20161015_750-50_EmptySlot_'
# nSubFiles = 116 *4 # NB actual files, not filenums. So 4x num raw files
# cubeInfoFile = 'cubeinfoFeb2017.idlvar'



saveFilePref = 'allSummedImsCube_20170227_maxmeth'
saveFilePref = 'WIP'

# Assumes HWP angles go as 0, 22.5, 45, 67.5
startFileNum = 0
# nSubFiles = 164 *4 # NB actual files, not filenums. So 4x num raw files
# nSubFiles = 172 *4 # NB actual files, not filenums. So 4x num raw files

showAllImages = False
showPlots = True
saveAllSummedIms = False

centroidThresh = 5000 #Values below this clipped to zero for centroiding
#centroidThresh = 2000

pixScale = 6.5 # Just for scale bar

######################################################################################################
# Settings for making the differential images

twoCamMode = True # If true, just reflect chan2 horizontally, not rotate.

#bsRot = 74.74 #Differential rotation of beamsplitter, in degrees
bsRot = 89.93 #Why is it not 74.74, like in make_tmpl????????


doMask = False
maskRad = 10 #Radius of central mask (px)
maskVal = 1  #Value to set central masked region to

flatFileC1 = ''
flatFileC2 = ''

clim=[-5, 5]
clim=[-100, 100]

dataPath = './'
dataPath = "/Users/bnorris/DontBackup/simplePDIdata/"
dataPath = "/Volumes/RedDisk/VAMPIRES_WorkingData/simplePDIdata/"
#dataPath = '/Volumes/BN_DATA_OSX/VAMPIRES_WorkingData/MWC758/'
# dataPath = '/Users/bnorris/DataAnalysis/VAMPIRES_DataAnalysis/SimplePDI/'
# dataPath = '/Users/bnorris/DontBackup/simplePDIdata/'

# loadFilename = 'allSummedImsCube_176files__cube_HD169142_01_20160721_675-50_EmptySlot_.npy'
# loadFilename = 'allSummedImsCube2_cube_HD169142_01_20160721_675-50_EmptySlot_CentroidLLim500.npy'
# loadFilename = 'allSummedImsCube_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
#loadFilename = 'allSummedImsCube_cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_.npy'
#loadFilename = 'allSummedImsCube_cube_HD19820_01_20160918_750-50_EmptySlot_.npz'
# loadFilename = 'allSummedImsCubeB_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
# loadFilename = 'allSummedImsCube_maxMeth_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
loadFilename = 'allSummedImsCube_centMeth_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
# loadFilename = 'allSummedImsCube_cenMeth_cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_.npz'
# loadFilename = 'allSummedImsCube_maxMeth_cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_.npz'
#loadFilename = dataPath+loadFilename
#loadFilename = 'allSummedImsCube_maxMeth_cube_AUMic20161015_01_Combined_20161015_750-50_EmptySlot_.npz'
loadFilename = 'allSummedImsCube_20170227cube_ABAur_01-02Combined_20161015_750-50_EmptySlot_.npz'
# loadFilename ='allSummedImsCube_20170227_maxmethcube_ABAur_01-02Combined_20161015_750-50_EmptySlot_.npz'
# loadFilename= 'allSummedImsCube_201702ABAur_01-03Combined_20170116_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_201701_332acqs_maxmethABAur_01-03Combined_20170116_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_201701_332acqs_COMmethABAur_01-03Combined_20170116_750-50_EmptySlot_0.npz'
loadFilename = 'allSummedImsCube_ABAur201612_52acqs_COMmethABAur_02_20161215_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_proc201709_COMmeth_ABAur_01_20170313_750-50_EmptySlot_0.npz'

## loadFilename = 'allSummedImsCube_201709_COMmeth_omiCet_01_20170912_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_201710_68acqs_COMmeth_subtDark2_omiCet_01_20170912_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_201710_68acqs_COMmeth_noBGsubt_omiCet_01_20170912_750-50_EmptySlot_0.npz'

# loadFilename = 'allSummedImsCube_201709_COMmeth_PREFLIP_omiCet_01_20170912_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_201711_COMmeth_subtdark_63Cet_01_20170912_750-50_EmptySlot_0.npz'
###loadFilename = 'allSummedImsCube_201704_128acqs_COMmethHD141569_02Combined_20170418_750-50_EmptySlot_0.npz'

doDerotate = True

# Only use a subset of images? This is in *HWP set* numbering
# imRange = [5, 41]
#imRange = [0, 29]

# Reload cube if allSummedIms already exists?
#reload = False


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
    def __init__(self):
        self.allSummedImsOrig = None
        self.allPAsOrig = None
        self.allStokesIms = None
        self.allPolzstateIms = None


    def loadCube(self, dataPath, loadFilename, twoCamMode = False):
        npzfile = np.load(dataPath+loadFilename)
        self.allSummedImsOrig = npzfile['arr_0']
        self.allPAsOrig = npzfile['arr_1']
        del (npzfile)
        print self.allPAsOrig

        if twoCamMode:
            print 'loadCube in twoCamMode'
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

        print 'Cube loaded.'


    def processRawFiles(self, dataPath, filePref, cubeInfoFile, nSubFiles, startFileNum = 0,
                          fileExtn = '.fits', saveFilePref = 'allSummedImsCube_', showImage = False,
                          centroidThresh=5000, saveCube = True, method = 'com', IDLPreCubed = False,
                          twoCamMode = False):

        if not IDLPreCubed and not twoCamMode:
            print('ERROR: Cannot do single-cam data without IDL pre-cubing.')
            return

        if IDLPreCubed:
            # Case when files are sub-cubes made from IDL pipeline
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
            # Note on PAs: As of Oct2017, the PA keyword in vampires data is just the IMR.PAD value, i.e. without the
            # pupil offset (IMR.PAP). To make the PA agree with iObserve values, it *seems* it is:
            # realPA = pad - pap - 180 for Northern targets, or
            # realPA = pad - pap + 180 for Southern targets, or
            cubeInfoFileString = dataPath + cubeInfoFile
            self.cubeInfo = self.readCubeInfo(cubeInfoFileString)
            pas = self.cubeInfo.pa
            # Indexes are [:, :, hwpSet, HWP, Channel, LCVRstate]
            # So each set of 4 VAMPIRES on-sky files corresponds to 1 hwpSet.

        else:
            # Case when files are raw camera files
            curFilenumStr = '%d' % startFileNum
            curCamStr = '_cam1'
            curFilename = dataPath + filePref + curFilenumStr + curCamStr + fileExtn
            hdulist = fits.open(curFilename)
            curHDU = hdulist[0]
            curCube = np.transpose(curHDU.data)
            pa = hdulist[1].header['PAD'] #TODO - Do proper PA rather than just PAD value
            hdulist.close()
            nFrms = curCube.shape[2]
            dim = curCube.shape[0]


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
                    print 'Reading file %s' % curFilename
                    print 'HWPSet %d -- HWP State %d, Channel %d, LCVR %d' % \
                          (curHWPSet, hwpState, curBSChan, curLCVRState)
                    hdulist = fits.open(curFilename)
                    curHDU = hdulist[0]
                    curCube = np.transpose(curHDU.data)
                    hdulist.close()

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
                imRange = [], showPlots = True, maskfile = [], calMat_PreRot = [],
                calMat_PostRot = [], twoCamMode = False, doFlats = False, Ibias = 0):
        if Ibias > 0:
            print('WARNING: Stokes I bias applied')
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
                                curImIn = np.fliplr(curImIn)
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
            curIm[:, :, 0] = I + Ibias # Ibias is a kludge to increase visibility of some features in fractional polz map
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



    def partialDiffcal(self, calmode, calMat_PreRot = []):
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
        print 'Done.'


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
                     radPower = 0, cmap = None, coeffs = None, norm=True, crop=None,
                     showColorbar=True):
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

        bLen = 100./pixScale
        plt.plot((10, 10+bLen), (10, 10), 'w')
        plt.text(10, 17, '~100 mas', color='w')
        # bLen = 100./pixScale
        # plt.plot((5, 5+bLen), (5, 5), 'w')
        # plt.text(5, 7, '~100 mas', color='w')

        plt.imshow(im, clim=clim, interpolation='Nearest', cmap=cmap)
        if showColorbar:
            plt.colorbar()

        print 'Image min: %f, max: %f' % (im.min(), im.max())
        return im




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




calMat = np.asarray([[1., 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])

# Includes no 90 deg rotation
calMat = np.asarray([[1.0653, 0.0457, 0.0084, 0],
                     [0.0228, 0.7153, -0.0355, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])
# No 90 deg, but included 3rd row before inv
calMat = np.asarray([[1.0655, 0.0458, 0.0060, 0],
                     [0.0220, 0.7147, -0.0254, 0],
                     [0.0233, 0.0163, 0.7139, 0],
                     [0, 0, 0, 1]])

# # Includes 90 deg rotation
# calMat = np.asarray([[1.0624, -0.0455, 0.0129, 0],
#                      [-0.0228, -0.7133, 0.0354, 0],
#                      [0, 0, 1, 0],
#                      [0, 0, 0, 1]])


#
p = pdiImages()

# p.processRawFiles(dataPath, filePref, cubeInfoFile, nSubFiles, startFileNum=startFileNum,
#                   saveFilePref = saveFilePref, centroidThresh = centroidThresh, showImage=False,
#                   method='max', twoCamMode = False, IDLPreCubed = True)
# p.loadCube(dataPath, loadFilename)
p.loadCube(dataPath, loadFilename, twoCamMode = twoCamMode)

# p.makeDiffIms(imRange=imRange, showPlots=False, twoCamMode=twoCamMode, deRotate=doDerotate)
p.makeDiffIms(showPlots=False, twoCamMode=twoCamMode, deRotate=doDerotate, Ibias=0)
#p.makeDiffIms(imRange=imRange, showPlots=False, maskfile='manualMasksNoc2016.npz')
#p.makeDiffIms(showPlots=False)
# mat1 = np.identity(4)
# mat2 = np.identity(4)
# p.makeDiffIms(imRange=imRange, showPlots=False, calMat_PreRot=mat1, calMat_PostRot=mat2)

# p.makeDiffIms(imRange=imRange, showPlots=False, twoCamMode=twoCamMode, deRotate=doDerotate, doFlats=True)

I = p.StokesImsSummed[:,:, 0]
Q = p.StokesImsSummed[:,:, 1]
U = p.StokesImsSummed[:,:, 2]
pQ = Q / I
pU = U / I
p2_frac = np.sqrt(pQ**2 + pU**2)
p_frac = np.sqrt(Q**2 + U**2) / I
P = np.sqrt(Q**2 + U**2)

#p.plotSingleIm(imType='p', clim=[-0.05, 0.4])
#p.plotSingleIm(imType='p', cmap='gray', crop=0.5)

# im = p2_frac
# im = np.clip(im,-0.1,0.4)
# # im = U
# #im = p.polzRatioImsSummed
# curHDU = fits.PrimaryHDU()
# curHDU.data = im
# curHDU.writeto('out.fits', clobber=True)

curHDU = fits.PrimaryHDU()
curHDU.data = np.transpose(p.StokesImsSummed)
curHDU.writeto('out.fits', clobber=True)


# plt.imshow(im)
# plt.imshow(im, cmap='gray')
p.plotSingleIm(imType='p',crop=0.5, cmap='viridis')
plt.pause(0.001)


# # Show polarisation ratio images
# im_f = p.polzRatioImsSummed
# plt.imshow(im_f)
# plt.imshow(im_f[48:208, 48:208])
# #plt.pause(0.001)

# fits.writeto('out.fits', data, header)

"""
dataPath = '/import/silo4/snert/VAMPIRESData_201609/Analysis/MWC758/'
filePref = 'cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_'
nSubFiles = 172 *4 # NB actual files, not filenums. So 4x num raw files

import simplePDIFull_01 as s
cubeInfoFile = 'cubeinfoOct2016.idlvar'
p = s.pdiImages()
p.processRawFiles(dataPath, filePref, cubeInfoFile, nSubFiles,saveFilePref =
    'allSummedImsCube_maxMeth_', method='max')
"""

print 'Done.'



"""
from rebin import rebin
rebinSize = 4

pQ=p.plotSingleIm(imType='pQ',crop=0.5, cmap='viridis')
pU=p.plotSingleIm(imType='pU',crop=0.5, cmap='viridis')

pQ_r = rebin(pQ, rebinSize)
pU_r = rebin(pU, rebinSize)
p_r = np.sqrt(pQ_r**2 + pU_r**2)
psi_r = 0.5*np.arctan(pU_r/pQ_r)
psi = 0.5*np.arctan(pU/pQ)

psi_r = -psi_r + np.pi/2 #0 should be up, not right...?

plt.clf()
plt.imshow(p_r)
# Use quiver but give the angles explicitly, and magnitudes by setting one input to fractional P
# plt.quiver(p_r, 0, angles=psi_r/np.pi*180, pivot='mid', headwidth=0, scale=0.4)
plt.quiver(p_r, 0, angles=psi_r/np.pi*180, pivot='mid', headwidth=0)
"""



