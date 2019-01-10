# Just align and stack the x% best images in a cube

path = '/Users/bnorris/DataAnalysis/201703_VAMPIRESStuff/'
filename = 'HD91312_02_20170313_750-50_EmptySlot_02_cam1.fits'
#filename='ABAur_01_20170313_750-50_EmptySlot_083_cam1.fits'
filename = 'procyon_01_Halpha_EmptySlot_0_cam1.fits'
filename = 'procyon_01zapoff_Halpha_EmptySlot_0_cam1.fits'

path = '/Users/bnorris/DontBackup/simplePDIdata/omiCet_01_20170912/'
path = '/Users/bnorris/data/simplePDIdata/omiCet_01_20170912/'
filename = 'omiCet_01_20170912_750-50_EmptySlot_00_cam1.fits'
# filename = '63Cet_01_20170912_750-50_EmptySlot_015_cam1.fits'

path = '/Users/bnorris/DontBackup/vampdata_Oct2018/'
filename = 'mwc758_01_20181023_Open_EmptySlot_020_cam1.fits'
filename = 'mwc758_02_20181023_Open_EmptySlot_031_cam1.fits'

title = 'With Zap - 16s Aligned & Summed'
title = 'Without Zap - 16s (best 20%) Aligned & Summed'
title = ''
cmap = 'viridis'
# cmap = 'jet'
# cmap = 'hot'

# Percentile cuts to use. List multiple values to do several runs
percents = (0, 20, 40, 60, 80, 90, 99)
percents = [90]
# percents = [80]

# Middle of a 256 cube
x1 = 96
x2 = 160
y1 = 96
y2 = 160

# # Cent 115,136
# x1 = 83
# x2 = 147
# y1 = 104
# y2 = 168

# # Middle of a 128 cube
# x1 = 32
# x2 = 96
# y1 = 32
# y2 = 96

# # Cent 54, 72
# x1 = 22
# x2 = 86
# y1 = 40
# y2 = 104

alignMethod = 'COM' # 'None', 'max' or 'COM'
boxRad = 4 # For COM method
makeMovie = False
saveCube = True


import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt
from copy import copy
import polzimtools as plz
import time
import matplotlib.animation as manimation
plt.ion()


cube = fits.getdata((path+filename))
cube = cube[1:-1, :, :]
nf = len(cube)
sz = cube.shape[1]

# Get the peak flux in each cube
pkfluxes = np.amax(cube,axis=(1,2))
#sortInds = np.argsort(pkfluxes)

stackedIms = []
# Make top x% images
for perc in percents:
    print(perc)
    cut = np.percentile(pkfluxes, perc)
    goodInds = np.where(pkfluxes >= cut)
    curCube = cube[goodInds[0], :, :]
    curCubeShifted = np.zeros_like(curCube)
    for c in range(curCube.shape[0]):
        curIm = curCube[c, :, :]
        imCent = sz / 2

        if alignMethod == 'max':
            # Max pixel method
            curCent = np.where(curIm == curIm.max())
        if alignMethod == 'COM':
            # Centre of mass method
            roughCent = np.where(curIm == curIm.max())
            mask = np.zeros_like(curIm)
            mask[roughCent[0][0]-boxRad:roughCent[0][0]+boxRad,
                roughCent[1][0] - boxRad:roughCent[1][0] + boxRad] = 1
            curImMasked = curIm * mask
            curCent = ndimage.center_of_mass(curImMasked)
            # plt.imshow(curImMasked)
            # plt.pause(0.1)

        if alignMethod == 'None':
            curCubeShifted[c, :, :] = curCube[c, :, :]
        else:
            rowOffset = imCent - curCent[0]
            colOffset = imCent - curCent[1]
            curCubeShifted[c, :, :] = ndimage.interpolation.shift(curCube[c, :, :],
              [rowOffset, colOffset])

    stackedIm = np.mean(curCubeShifted, axis=0)
    stackedIm = stackedIm[y1:y2, x1:x2]
    stackedIms.append(stackedIm)
    plt.imshow(stackedIm, interpolation = 'nearest', cmap=cmap)
    plt.title(title)
    plt.draw()
    plt.show()
    # plt.pause(1)


if saveCube:
    a = curCubeShifted[:, y1:y2, x1:x2]
    stackedImsArr = np.asarray(stackedIms)
    fits.writeto('cubeOut.fits', a, clobber=True)
    fits.writeto('stackedImOut.fits', stackedIm, clobber=True)
    fits.writeto('stackedImsOut.fits', stackedImsArr, clobber=True)

if makeMovie:
    # I can't get this to work:
    # FFMpegWriter = manimation.writers['ffmpeg']
    # metadata = dict(title='SCExAOVAMPIRES PSF Demo', artist='',
    #                 comment='')
    # writer = FFMpegWriter(fps=25, metadata=metadata)
    # fig = plt.figure(2)
    # l, = plt.imshow(stackedIms[0], interpolation='nearest')

    nFrms = len(curCubeShifted)

    for i in range(nFrms):
        plt.clf()
        plt.imshow(curCubeShifted[i,y1:y2, x1:x2], interpolation = 'nearest', cmap='viridis')
        # plt.imshow(np.log10(curCubeShifted[i, y1:y2, x1:x2]), interpolation='nearest', cmap=cmap)
        filename = './movie_out/imsave5_' + str(i) + '.png'
        plt.savefig(filename, bbox_inches='tight')


