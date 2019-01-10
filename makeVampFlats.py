"""
This will assemble flats from FLC acqusition sequences, for the 2-camera setup.
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

saveCube = True
saveFilePref = 'flats_noPolz_107-160_badPixCorr-smoothed20_'
saveFilePref = 'flats_noPolz_107-160_fullDarkSub_badPixCorr-smoothed20_'

# Assumes HWP angles go as 0, 22.5, 45, 67.5
startFileNum = 0

# What to do about bg subtraction?
# 0 - don't subtract anything
# 1 - subtract specified bias values, as per bgVals tuple
# 2 - subtract actual darkframes from file specified in darkFilename
# 3 - use a sky annulus, as per skyAnnInner and Outer values. (Careful - halo is big!)-NOT IMPL.
bgMode = 2
bgVals = (176.5, 178.9) # (chan1, chan2)
darkFilename = 'summedDarks_darks_20170118_512_em1000_18000us_20170119_750-50_Mirror_0.npz'

# Remove bad cam2 pixels -
# This is to get rid of the problem pixels in the faulty early 2017 camera 2
# All settings are in the actual code...
removeBadCam2 = True

# Gaussian filter result - sigma~10 removes most pixel-scale noise. 0 for none.
gaussFilter = True
gaussFilterSig = 20

showPlots = True


# Skip reading raw files and instead read previously saved cube
readPrevCube = False
cubePath = '/Users/bnorris/DontBackup/simplePDIdata/'
cubeName = 'flats_noPolz_107-160_domeflat_nopolz_05_20170117_750-50_EmptySlot_0.npz'



#####################################################################################################

if bgMode == 2:
    npzfile = np.load(darkFilename)
    darkFrames = npzfile['finalDarks']
    del (npzfile)

if readPrevCube:
    npzfile = np.load(cubePath + cubeName)
    finalFlats = npzfile['arr_0']

else:
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

                if bgMode == 2:
                    curCubeSummed = curCubeSummed - darkFrames[:, :, c]
                    darkMed = np.median(darkFrames[:, :, c])
                    print('Subtracting dark frame (with median %f)' % darkMed)

                print 'Max val in summed unshifted image: %f' % np.max(curCubeSummed.ravel())

                allSummedIms[:, :, curHWPSet, hwpState, c, l] = curCubeSummed


        # Increment HWP state
        hwpState = hwpState + 1
        if hwpState == 4:
            hwpState = 0
            curHWPSet = curHWPSet + 1


    # Sum the subcubes
    finalFlats = np.mean(allSummedIms, axis=2)

if removeBadCam2:

    # Good for original bgVals = (176.5, 178.9) subtraction
    badpix_llim = 240
    badpix_hlim = 320 #300

    # Seems to work ok for darkframe-subtracted
    badpix_llim = 240
    badpix_hlim = 300

    finalFlatsNobadpx = copy(finalFlats)
    for h in range(0, 4):
        for l in range(0, 2):
            im_c2 = copy(finalFlats[:, :, h, 1, l])
            mask_c2 = np.ones_like(im_c2)
            mask_c2[((im_c2 < badpix_llim) | (im_c2 > badpix_hlim))] = 0

            mask_c2[425:, :] = 1
            mask_c2[:, :80] = 1

            medImage = ndimage.median_filter(im_c2, size=8)
            medImage = ndimage.gaussian_filter(medImage, sigma=16)
            im_c2[mask_c2 == 0] = medImage[mask_c2 == 0]
            finalFlatsNobadpx[:, :, h, 1, l] = im_c2

            if showPlots:
                plt.clf()
                plt.subplot(1, 2, 1)
                # plt.imshow(finalFlatsNobadpx[:, :, h, 0, l], interpolation='nearest')
                plt.imshow(finalFlatsNobadpx[:, :, h, 0, l], interpolation='nearest', clim = (160, 330))
                plt.colorbar()
                plt.subplot(1, 2, 2)
                plt.imshow(finalFlatsNobadpx[:, :, h, 1, l], interpolation='nearest', clim = (160, 330))
                plt.colorbar()
                plt.pause(0.2)

    finalFlats = finalFlatsNobadpx
    plt.pause(1)

# if bgMode == 2:
#     for h in range(0, 4):
#         for l in range(0, 2):
#             for c in range(0, 2):
#                 finalFlats[:, :, h, c, l] = finalFlats[:, :, h, c, l] - darkFrames[:, :, c]
#                 darkMed = np.median(darkFrames[:, :, c])
#                 print('Subtracting dark frame (with median %f)' % darkMed)

if gaussFilter:
    for h in range(0, 4):
        for l in range(0, 2):
            for c in range(0, 2):
                finalFlats[:, :, h, c, l] = ndimage.gaussian_filter(finalFlats[:, :, h, c, l],
                                                                    sigma = gaussFilterSig)
            if showPlots:
                plt.clf()
                plt.subplot(1, 2, 1)
                plt.imshow(finalFlats[:, :, h, 0, l], clim=(160, 330))
                plt.colorbar()
                plt.subplot(1, 2, 2)
                plt.imshow(finalFlats[:, :, h, 1, l], clim=(160, 330))
                plt.colorbar()
                plt.pause(0.2)


plt.clf()
plt.subplot(1, 2, 1)
plt.imshow(finalFlats[:, :, 0, 0, 0], interpolation='nearest')
plt.colorbar()
plt.subplot(1, 2, 2)
plt.imshow(finalFlats[:, :, 0, 1, 0], interpolation='nearest')
plt.colorbar()

if saveCube:
    saveFilename = saveFilePref + filePref
    np.savez(saveFilename, finalFlats)








