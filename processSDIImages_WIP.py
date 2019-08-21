import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt
from copy import copy

fileExtn = '.fits'
dataPath = "/Users/bnorris/DontBackup/simplePDIdata/Ha/"
dataPath = "/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201709/20170911/"

filePref = 'omiCet_01_HaDifferential-'
fileSuf = '_Open_EmptySlot_'
startFileNum = 0
nSets = 8 #1 # 1 Set is two states x two cams

saveFilePref = 'allSummedImsCube_WIP_'
centroidThresh = 2000
showAllImages = False
showPlots = True
saveAllSummedIms = False
method = 'com' # 'com' or 'max'
saveCube = False #True
skipFirstFrame = True

# Filter arrangement, listed as [cam1, cam2]
state1 = ('cont', 'Ha')
state2 = ('Ha', 'cont')


# e.g. filename: omiCet_01_HaDifferential-state1_Open_EmptySlot_0_cam1.fits
#################################################################################################

# Read a FITS file to get sizes
curFilenumStr = '%d' % startFileNum
curCamStr = '_cam1'
curStateStr = 'state1'
curFilename = dataPath + filePref + curStateStr + fileSuf + curFilenumStr + curCamStr + fileExtn
hdulist = fits.open(curFilename)
curHDU = hdulist[0]
curCube = np.transpose(curHDU.data)
nFrms = curCube.shape[2]
if skipFirstFrame:
    nFrms = nFrms - 1
dim = curCube.shape[0]


# Indexes are [:, :, Set, State, Channel (camera)]
allSummedIms = np.zeros([dim, dim, nSets, 2, 2])
pas = [] # List of PAs from each file read
maxvals = [] # For testing

for f in range(0, nSets):
    curFileNum = f
    curFilenumStr = '%d' % curFileNum

    for c in range(0, 2):
        curChan = c
        if curChan == 0:
            curCamStr = '_cam1'
        else:
            curCamStr = '_cam2'

        for s in range(0, 2):
            curState = s
            if curState == 0:
                curStateStr = 'state1'
            else:
                curStateStr = 'state2'

            curFilename = dataPath + filePref + curStateStr + fileSuf + curFilenumStr \
                    + curCamStr + fileExtn
            print('Reading file %s' % curFilename)
            hdulist = fits.open(curFilename)
            curHDU = hdulist[0]
            curCube = np.transpose(curHDU.data)
            if skipFirstFrame:
                curCube = curCube[:, :, 1:]

            # Get PA from header (use PAD, and assume instrument offset added later)
            # Only do for each new state (since same for 2 simul cams)
            if c == 0:
                # pa = hdulist[1].header['PAD']
                pa = 0
                pas.append(pa)

            hdulist.close()

            if c == 1:
                curCube = np.fliplr(curCube)

            curCubeSummed = np.mean(curCube, axis=2)
            print('Max val in summed UNshifted image: %f' % np.max(curCubeSummed.ravel()))

            # Make cube of shifted images
            curCubeShifted = np.zeros([dim, dim, nFrms])

            # Centre the image
            for i in range(0, nFrms):

                llim = centroidThresh
                curIm = copy(curCube[:, :, i])

                if method is 'com':
                    curIm[np.where(curIm <= llim)] = 0
                    if np.sum(curIm.ravel()) <= 0:
                        print('WARNING: Frame found with no vals above COM threshold')
                    curCent = ndimage.center_of_mass(curIm)
                elif method is 'max':
                    curCent = np.where(curIm == curIm.max())
                else:
                    print('Error: Unknown centering method specified')

                if showAllImages:
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

            allSummedIms[:, :, f, s, c] = curCubeShiftedSummed
            maxvals.append(np.max(curCubeShiftedSummed.ravel()))
pas = np.asarray(pas)

if saveCube:
    # Save the big cube
    saveFilename = saveFilePref + filePref
    np.savez(saveFilename, allSummedIms, pas)




##### Now go through and make a list of difference images
# SingleDiffImages arranged such that each entry is Ha - cont
allDoubleDiffImages = np.zeros(([dim, dim, nSets]))
allSingleDiffImages = np.zeros(([dim, dim, nSets, 2]))
for k in range(0, nSets):
    allSingleDiffImages[:, :, k, 0] = allSummedIms[:, :, k, 0, 1] - allSummedIms[:, :, k, 0, 0]
    # allSingleDiffImages[:, :, k, 1] = allSummedIms[:, :, k, 1, 0] - allSummedIms[:, :, k, 1, 1]
    allSingleDiffImages[:, :, k, 1] = allSummedIms[:, :, k, 1, 1] - allSummedIms[:, :, k, 1, 0]
    allDoubleDiffImages[:, :, k] = (allSummedIms[:, :, k, 0, 1] - allSummedIms[:, :, k, 0, 0]) \
        - (allSummedIms[:, :, k, 1, 1] - allSummedIms[:, :, k, 1, 0])

summedDoubDiffIm = np.mean(allDoubleDiffImages, axis=2)
summedSingleDiffIm1 = np.mean(allSingleDiffImages[:, :, :, 0], axis=2)
summedSingleDiffIm2 = np.mean(allSingleDiffImages[:, :, :, 1], axis=2)



##### Try differnces of normalised images
ab = allSummedIms[:, :, k, 0, 1] / np.max(allSummedIms[:, :, k, 0, 1])
aa = allSummedIms[:, :, k, 0, 0] / np.max(allSummedIms[:, :, k, 0, 0])
bb = allSummedIms[:, :, k, 1, 1] / np.max(allSummedIms[:, :, k, 1, 1])
ba = allSummedIms[:, :, k, 1, 0] / np.max(allSummedIms[:, :, k, 1, 0])
dIm = (ab - aa) - (bb - ba)


##### Try ratio images
k = 0
#allDoubleRatioImages = np.zeros(([dim, dim, nSets]))
r1 = (allSummedIms[:, :, k, 0, 1] / allSummedIms[:, :, k, 0, 0])
r2 = (allSummedIms[:, :, k, 1, 1] / allSummedIms[:, :, k, 1, 0])
# r1 = (allSummedIms[:, :, k, 1, 0] / allSummedIms[:, :, k, 1, 1])
# r2 = (allSummedIms[:, :, k, 0, 0] / allSummedIms[:, :, k, 0, 1])
R = np.sqrt(r1/r2)
rIm = (R-1)/(R+1)



plt.imshow((dIm)[32:224, 32:224], cmap='viridis')
plt.imshow((rIm)[32:224, 32:224], cmap='viridis')

plt.imshow((R)[96:160, 96:160], cmap='viridis')
plt.imshow((allSummedIms[:, :, k, 0, 1])[96:160, 96:160], cmap='viridis')
