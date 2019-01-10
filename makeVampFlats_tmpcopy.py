"""
This will assemble flats from FLC acqusition seuqences, for the 2-camera setup.
"""


import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt
from copy import copy
import polzimtools as plz
import time

fileExtn = '.fits'
dataPath = "/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201701/20170117/"
filePref = 'domeflat_nopolz_05_20170117_750-50_EmptySlot_0'

# NB separate subfiles for each state, not original filenums.
# So for single camera (old BS), this is 4x num raw files (because 4 states per file)
# For dual-camera, this is 2x num raw files (because 2 states (FLC) per file)
#   And num raw files is 2x the maximum numbered file (since 2 camera files)
#   i.e. for files numbered 0,..,7 , nSubFiles = 8 * 2 * 2 = 32
nSubFiles = 24 *2*2
# nSubFiles = 4 *2*2

saveFilePref = 'flats_20171105_noPolz_107-160_24sets'

# Assumes HWP angles go as 0, 22.5, 45, 67.5
startFileNum = 0

# What to do about bg subtraction?
# 0 - don't subtract anything
# 1 - subtract specified bias values, as per bgVals tuple
# 2 - subtract actual darkframes from file specified in darkFilename - NOT IMPLEMENTED
# 3 - use a sky annulus, as per skyAnnInner and Outer values. (Careful - halo is big!)-NOT IMPL.
bgMode = 1
bgVals = (176.5, 178.9) # (chan1, chan2)





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

        # Now sort frames into the 2 FLC states and discard 1st 2 frames
        nSubFrms = nFrms / 2 - 1
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
            print 'HWPSet %d -- HWP State %d, Channel %d, LCVR %d' % \
                  (curHWPSet, hwpState, curChan, curFLCState)

            if l == 0:
                curCube = curCube_FLC1
            else:
                curCube = curCube_FLC2

            curCubeSummed = np.mean(curCube, axis=2)

            if bgMode == 1:
                curCubeSummed = curCubeSummed - bgVals[curChan]
                print('Dark value subtracted: %f' % bgVals[curChan])

            print 'Max val in summed unshifted image: %f' % np.max(curCubeSummed.ravel())

            allSummedIms[:, :, curHWPSet, hwpState, c, l] = curCubeSummed


    # Increment HWP state
    hwpState = hwpState + 1
    if hwpState == 4:
        hwpState = 0
        curHWPSet = curHWPSet + 1


# Sum the subcubes
finalFlats = np.mean(allSummedIms, axis=2)


plt.figure()
plt.imshow(finalFlats[:, :, 0, 1, 0])

# Save the big cube
saveFilename = saveFilePref + filePref
np.savez(saveFilename, allSummedIms)
















