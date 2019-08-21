import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt


dataPath = '/Users/bnorris/DontBackup/vampFlats/'

inFiles = ['domeflats_750_1s_em300_750-50_EmptySlot_0_cam1.fits', 'domeflats_750_1s_em300_750-50_EmptySlot_1_cam1.fits',
           'domeflats_750_1s_em300_750-50_EmptySlot_2_cam1.fits', 'domeflats_750_1s_em300_750-50_EmptySlot_3_cam1.fits']
outFile = 'domeflats_20170419_1s_em300_750_cam1.fits'

inFiles = ['domeflats_750_1s_em300_750-50_EmptySlot_0_cam2.fits', 'domeflats_750_1s_em300_750-50_EmptySlot_1_cam2.fits',
           'domeflats_750_1s_em300_750-50_EmptySlot_2_cam2.fits', 'domeflats_750_1s_em300_750-50_EmptySlot_3_cam2.fits']
outFile = 'domeflats_20170419_1s_em300_750_cam2.fits'




allIms = []
for file in inFiles:
    cube = fits.getdata(dataPath+file)
    cube = np.transpose(cube)
    im = np.mean(cube, axis=2)
    allIms.append(im)
    plt.clf()
    plt.imshow(im)
    plt.colorbar()
    plt.pause(1)

allImsArr = np.asarray(allIms)
finalIm = np.mean(allImsArr, axis=0)
plt.clf()
plt.imshow(finalIm)
plt.colorbar()

fits.writeto(dataPath+outFile, finalIm)