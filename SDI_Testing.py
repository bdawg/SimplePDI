import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt
from copy import copy

dataPath = "/Users/bnorris/DontBackup/simplePDIdata/Ha/"
filename = 'allSummedImsCube_omiCet_Combined01-04_HaDifferential.npz'

# Filter arrangement, listed as [cam1, cam2]
state1 = ('cont', 'Ha')
state2 = ('Ha', 'cont')


npzfile = np.load(dataPath+filename)
allSummedIms = npzfile['arr_0']
allPAs = npzfile['arr_1']

nSets = allSummedIms.shape[2]
dim = allSummedIms.shape[0]

allDoubleRatioImages = np.zeros(([dim, dim, nSets]))
allDoubleDiffNormImages = np.zeros(([dim, dim, nSets]))
allContImages = np.zeros(([dim, dim, nSets]))
allHaImages = np.zeros(([dim, dim, nSets]))

for k in range(0, nSets):
    r1 = (allSummedIms[:, :, k, 0, 1] / allSummedIms[:, :, k, 0, 0])
    r2 = (allSummedIms[:, :, k, 1, 1] / allSummedIms[:, :, k, 1, 0])
    # r1 = (allSummedIms[:, :, k, 1, 0] / allSummedIms[:, :, k, 1, 1])
    # r2 = (allSummedIms[:, :, k, 0, 0] / allSummedIms[:, :, k, 0, 1])
    R = np.sqrt(r1 / r2)
    rIm = (R - 1) / (R + 1)
    # allDoubleRatioImages[:, :, k] = rIm
    allDoubleRatioImages[:, :, k] = R

    contIm = (allSummedIms[:, :, k, 0, 0] + allSummedIms[:, :, k, 1, 1]) /2
    allContImages[:, :, k] = contIm
    HaIm = (allSummedIms[:, :, k, 0, 1] + allSummedIms[:, :, k, 1, 0]) /2
    allHaImages[:, :, k] = HaIm

    ab = allSummedIms[:, :, k, 0, 1] / np.max(allSummedIms[:, :, k, 0, 1])
    aa = allSummedIms[:, :, k, 0, 0] / np.max(allSummedIms[:, :, k, 0, 0])
    bb = allSummedIms[:, :, k, 1, 1] / np.max(allSummedIms[:, :, k, 1, 1])
    ba = allSummedIms[:, :, k, 1, 0] / np.max(allSummedIms[:, :, k, 1, 0])
    dIm = (ab - aa) - (bb - ba)
    allDoubleDiffNormImages[:, :, k] = dIm



summedDoubleRatioImages = np.mean(allDoubleRatioImages, axis=2)
R = summedDoubleRatioImages
dIm = np.mean(allDoubleDiffNormImages, axis=2)
summedCont = np.mean(allContImages, axis=2)
summedHa = np.mean(allHaImages, axis=2)
summedContN = summedCont / np.max(summedCont)

plt.imshow(summedCont[96:160, 96:160], cmap='viridis')
plt.imshow(R[96:160, 96:160], cmap='viridis')
plt.imshow((R-summedContN)[96:160, 96:160], cmap='viridis')

# fits.writeto('autosave_R.fits', R, clobber=True)
# fits.writeto('autosave_cont.fits', summedCont, clobber=True)