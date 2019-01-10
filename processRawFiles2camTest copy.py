import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt
from copy import copy
import polzimtools as plz
import time

fileExtn = '.fits'
dataPath = "/Volumes/GLINTDATA/VAMPIRESDataCopies/VAMPIRESData_201612/20161215/"
#filePref = 'ABAur_01_20170116_750-50_EmptySlot_0'
filePref = 'ABAur_02_20161215_750-50_EmptySlot_0'

# NB separate subfiles for each state, not original filenums.
# So for single camera (old BS), this is 4x num raw files (because 4 states per file)
# For dual-camera, this is 2x num raw files (because 2 states (FLC) per file)
#   And num raw files is 2x the maximum numbered file (since 2 camera files)
#   i.e. for files numberd 0,..,7 , nSubFiles = 8 * 2 * 2 = 32
nSubFiles = 52 *2*2 #68 *4

###cubeInfoFile = 'cubeinfoFeb2017.idlvar'
saveFilePref = 'allSummedImsCube_ABAur201612_52acqs_MAXmeth'

# Assumes HWP angles go as 0, 22.5, 45, 67.5
startFileNum = 0

centroidThresh=5000
showImage = False
method = 'max'
saveCube = True
showPlots = True


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

allSummedIms = np.zeros([dim, dim, nSubFiles / 16, 4, 2, 2])
allRawPAs = [] # List of PAs from each file read
# curFileNum = startFileNum
hwpState = 0  # This counts through 4 positions
curHWPSet = 0  # Increments for each new set of HWP posns

for f in range(0, nSubFiles / 4):
    curFileNum = f

    for c in range(0, 2):
        curChan = c

        # Generate current filename
        if curChan == 0:
            curCamStr = '_cam1'
        else:
            curCamStr = '_cam2'

        curFilenumStr = '%d' % curFileNum
        curFilename = dataPath + filePref + curFilenumStr + curCamStr + fileExtn

        print 'Reading file %s' % curFilename
        # print 'HWPSet %d -- HWP State %d, Channel %d, LCVR %d' % \
        #       (curHWPSet, hwpState, curChan, curFLCState)
        hdulist = fits.open(curFilename)
        curHDU = hdulist[0]
        curSuperCube = np.transpose(curHDU.data) # 'Super' since both FLC states

        # Get PA from header (use PAD, and assume instrument offset added later)
        pa = hdulist[1].header['PAD']
        allRawPAs.append(pa)

        # Now sort frames into the 2 FLC states and discard 1st 2 frames
        nSubFrms = nFrms / 2 - 1
        curCube_FLC1 = np.zeros((dim, dim, nSubFrms))
        curCube_FLC2 = np.zeros((dim, dim, nSubFrms))
        curFLC = 1
        count = 0
        for k in range(2, nFrms):
            if curFLC == 1:
                curCube_FLC1[:, :, count] = curSuperCube[:, :, k]
                # print " "
                # print "curCube_FLC1:"
                # print k
                # print count
                # print curFLC

            else:
                curCube_FLC2[:, :, count] = curSuperCube[:, :, k]
                # print " "
                # print "curCube_FLC2:"
                # print k
                # print count
                # print curFLC
                count = count + 1

            if curFLC == 1:
                curFLC = 2
            else:
                curFLC = 1


        for l in range(0, 2):
            curFLCState = l
            print 'HWPSet %d -- HWP State %d, Channel %d, LCVR %d' % \
                  (curHWPSet, hwpState, curChan, curFLCState)

            if l == 0:
                curCube = curCube_FLC1
            else:
                curCube = curCube_FLC2

            curCubeSummed = np.mean(curCube, axis=2)
            print 'Max val in summed UNshifted image: %f' % np.max(curCubeSummed.ravel())

            # Make cube of shifted images
            curCubeShifted = np.zeros_like(curCube)
            # Centre the image
            for i in range(0, nSubFrms):
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

# This is to maintain compatability with the single-camera data
pas = np.repeat(allRawPAs, 2)

if saveCube:
    # Save the big cube
    saveFilename = saveFilePref + filePref
    np.savez(saveFilename, allSummedIms, pas)





