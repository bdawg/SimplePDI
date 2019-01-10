import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt
from copy import copy
import polzimtools as plz
import time

fileExtn = '.fits'

dataPath = "/Volumes/GLINTDATA/VAMPIRESDataCopies/VAMPIRESData_201701/20170116/linked/"
#filePref = 'ABAur_01_20170116_750-50_EmptySlot_0'
filePref = 'ABAur_01-03Combined_20170116_750-50_EmptySlot_0'

# NB separate subfiles for each state, not original filenums.
# So for single camera (old BS), this is 4x num raw files (because 4 states per file)
# For dual-camera, this is 2x num raw files (because 2 states (FLC) per file)
#   And num raw files is 2x the maximum numbered file (since 2 camera files)
#   i.e. for files numberd 0,..,7 , nSubFiles = 8 * 2 * 2 = 32
nSubFiles = 332 *2*2 #68 *4


dataPath = '/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201704/20170418/linkCopied/'
filePref = 'HD141569_02Combined_20170418_750-50_EmptySlot_0'
nSubFiles = 128 *2*2

saveFilePref = 'allSummedImsCube_201701_332acqs_COMmeth'
saveFilePref = 'allSummedImsCube_201704_128acqs_COMmeth'

dataPath = "/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201709/20170912/"
dataPath = '/Users/bnorris/DontBackup/simplePDIdata/omiCet_01_20170912/'
filePref = 'omiCet_01_20170912_750-50_EmptySlot_0'
nSubFiles = 68 *2*2
saveFilePref = 'allSummedImsCube_201710_68acqs_COMmeth_bothsubtchan1vals_'

dataPath = "/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201709/20170912/"
filePref = '63Cet_01_20170912_750-50_EmptySlot_0'
nSubFiles = 16 *2*2
saveFilePref = 'allSummedImsCube_201711_MAXmeth_subtdark_'


dataPath = "/Users/bnorris/DontBackup/vampdata_Oct2018/"
filePref = 'ABAur_01_20181018_750-50_EmptySlot_0'
nSubFiles = 40 *2*2
saveFilePref = 'allSummedImsCube_MAXmeth_darksubt_'


# Assumes HWP angles go as 0, 22.5, 45, 67.5
startFileNum = 0

centroidThresh=5000
showImage = False
method = 'max'
saveCube = True
showPlots = True


# What to do about bg subtraction?
# 0 - don't subtract anything
# 1 - subtract specified bias values, as per bgVals tuple
# 2 - subtract actual darkframes from file specified in darkFilename
# 3 - use a sky annulus, as per skyAnnInner and Outer values. (Careful - halo is big!)
bgMode = 2
bgVals = (176.5, 178.9) # (chan1, chan2)
darkFilename = 'summedDarks_darks_10ms_em200_20181017_750-50_Mirror_0.npz'
# Be careful with outer radius to avoid satellite spots!
skyAnnInner = 90
skyAnnOuter = 120



if bgMode == 2:
    npzfile = np.load(darkFilename)
    darkFrames = npzfile['finalDarks']
    del (npzfile)

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

allSummedIms = np.zeros([dim, dim, nSubFiles // 16, 4, 2, 2])
allRawPAs = [] # List of PAs from each file read
# curFileNum = startFileNum
hwpState = 0  # This counts through 4 positions
curHWPSet = 0  # Increments for each new set of HWP posns

for f in range(0, nSubFiles // 4):
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

        print('Reading file %s' % curFilename)
        # print 'HWPSet %d -- HWP State %d, Channel %d, LCVR %d' % \
        #       (curHWPSet, hwpState, curChan, curFLCState)
        hdulist = fits.open(curFilename)
        curHDU = hdulist[0]
        curSuperCube = np.transpose(curHDU.data) # 'Super' since both FLC states

        # Get PA from header (use PAD, and assume instrument offset added later)
        pa = hdulist[1].header['PAD']
        allRawPAs.append(pa)

        # Now sort frames into the 2 FLC states and discard 1st 2 frames
        nSubFrms = nFrms // 2 - 1
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
            print('HWPSet %d -- HWP State %d, Channel %d, LCVR %d' % \
                  (curHWPSet, hwpState, curChan, curFLCState))

            if l == 0:
                curCube = curCube_FLC1
            else:
                curCube = curCube_FLC2

            curCubeSummed = np.mean(curCube, axis=2)
            print('Max val in summed UNshifted image: %f' % np.max(curCubeSummed.ravel()))

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
                    print('Error: Unknown centering method specified')
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

                if bgMode == 2: #Darkframe subtraction must happen before shifting
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
                curHDU.writeto('autosave.fits', clobber=True)
            except:
                print("Warning: Couldn't autosave")

            if showPlots:
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

            if bgMode == 1:
                curCubeShiftedSummed = curCubeShiftedSummed - bgVals[curChan]
                print('Dark value subtracted: %f' % bgVals[curChan])

            if bgMode == 2:
                print('Dark frame subtracted, with median %f' % darkMed)

            if bgMode == 3:
                px = np.arange(dim) - dim / 2
                xv, yv = np.meshgrid(px, px)
                rpix = np.sqrt(xv ** 2 + yv ** 2)
                skyInds = np.where(np.logical_and(rpix >= skyAnnInner, rpix <= skyAnnOuter))
                skyVals = curCubeShiftedSummed[skyInds]
                skyVal = np.median(skyVals)
                print('Sky value subtracted: %f' % skyVal)
                curCubeShiftedSummed = curCubeShiftedSummed - skyVal

            allSummedIms[:, :, curHWPSet, hwpState, c, l] = curCubeShiftedSummed

    #         break
    #     break
    # break



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




"""
#Test code for sky bg:
imForSky = curCubeShiftedSummed

# Set inner and outer radii for sky annulus. Refer to a log stretch of image.
# e.g. plt.imshow(np.log10(curCubeShiftedSummed), clim=[2, 3])
# Be careful with outer radius to avoid satellite spots!
skyAnnInner = 90
skyAnnOuter = 120

px = np.arange(dim)-dim/2
xv, yv = np.meshgrid(px, px)
rpix = np.sqrt(xv**2 + yv**2)
skyInds = np.where(np.logical_and(rpix >= skyAnnInner, rpix <= skyAnnOuter))
skyVals = imForSky[skyInds]

tmpIm = np.zeros([dim, dim])
tmpIm[skyInds] = imForSky[skyInds]
plt.clf()
plt.imshow(tmpIm)
"""
