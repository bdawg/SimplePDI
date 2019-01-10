from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg.linalg import multi_dot as md
from time import clock as cl
import timeit
import numba
from scipy import ndimage

import polzimtools as plz


# Open a Stokes image:
infile = 'toyDisk.fits'

hdulist = fits.open(infile)
curHDU = hdulist[0]
inIm = np.transpose(curHDU.data, (1, 2, 0))
inIm = inIm.astype('f8')

# t = cl()
# allIms = plz.StokesToPol(inIm)
# print cl()-t



# Try rotating polz of images...
paRange = np.arange(-45, 45, 10)
allRotatedImages = []
allRotatedPolIms = []
imsz = inIm.shape[0]

for angle in paRange:
    print "Doing angle %f" % angle
    curIm = plz.rotImPolz(inIm, angle)
    allRotatedImages.append(curIm)

    curPolIm = plz.stokesToPol(curIm)
    allRotatedPolIms.append(curPolIm)

    #plt.imshow(curIm[:,:,1])
    plt.imshow(curPolIm[0])
    plt.pause(0.001)

print 'Done.'


