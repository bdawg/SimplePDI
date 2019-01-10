# Extract the appropriate flats form flatfield files

# From the IDL pipeline:
# speckpos = [[146,356],$ ;Manually set position of speckle centers, in form
#             [368,154]]  ; [ [Xup, Yup],[Xdown,YDown]]


posn1 = [146, 356]
posn2 = [368, 154]
imsz = 256
zoomSz = 128 # Smaller than actual imsz, for viewing

dataPath = "/Users/bnorris/DontBackup/simplePDIdata/"
rawFlatfileName = 'domeflat_750_hwp0_750-50_EmptySlot_0.fits'
rawFlatfileName = 'domeflat_imr90_750-50_EmptySlot_0.fits'

dataPath = './'
rawFlatfileName= 'summedNPFlat_Aug2015_newBS_750_IMR90.fits'
rawFlatfileName = 'flatframe_750_Aug2015.fits'

smoothAmt = 20

makeManually = False


#For checking alignment if needed:
testimage1Name = 'ABAur/cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_00000_1_A.fits'
testimage2Name = 'ABAur/cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_00000_2_A.fits'
testRawImageName = 'ABAur_01_20160918_750-50_EmptySlot_00.fits'


import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt
from copy import copy
import polzimtools as plz
import time

# rawFlatCube = fits.getdata(dataPath+rawFlatfileName)
# rawFlat = rawFlatCube.mean(axis=0)

rawFlat = fits.getdata(dataPath+rawFlatfileName)



# # Check how noisy these flats are:
# f1=np.mean(rawFlatCube[0:14,:,:],axis=0)
# f2=np.mean(rawFlatCube[15:-1,:,:],axis=0)
# f1 = ndimage.filters.gaussian_filter(f1,smoothAmt)
# f2 = ndimage.filters.gaussian_filter(f2,smoothAmt)
# diff=f1-f2
# print np.std(diff)
# Answer for domeflats: very. Looks like photon noise (dark bits less difference).
#   Mean flat value is ~7000 and std of difference is ~350.
# If smoothed by 5 pixel gaussian, get std = 21 over *whole* raw image

# Try smoothing with a gaussian... ok over the regions of interest maybe?
flatSmoothed = ndimage.filters.gaussian_filter(rawFlat,smoothAmt)


# # Check subframe alignment...
# # Note that this seems to leave a residual between the extracted subframes and
# # the output of the IDl pipeline by ~1e-4 units. Since these are supposedly integer
# # pixel values this is hopefully negligible (some numerical thing).
# testCube1 = fits.getdata(dataPath+testimage1Name)
# testCube2 = fits.getdata(dataPath+testimage2Name)
# testRawCube = fits.getdata(dataPath+testRawImageName)
# testim1 = testCube1[0,:,:]
# testim2 = testCube2[0,:,:]
# testRawIm = testRawCube[1,:,:]
# w = imsz/2
# testExtr1 = testRawIm[posn1[1]-w:posn1[1]+w, posn1[0]-w:posn1[0]+w]
# testExtr2 = testRawIm[posn2[1]-w:posn2[1]+w, posn2[0]-w:posn2[0]+w]


# Extract subframe flats from averaged raw flat
w = imsz/2
subflat1 = rawFlat[posn1[1]-w:posn1[1]+w, posn1[0]-w:posn1[0]+w]
subflat2 = rawFlat[posn2[1]-w:posn2[1]+w, posn2[0]-w:posn2[0]+w]
print "subflat1 xrange = %d to %d, yrange = %d to %d" % (posn1[1]-w,
    posn1[1]+w, posn1[0]-w, posn1[0]+w)
print "subflat2 xrange = %d to %d, yrange = %d to %d" % (posn2[1]-w,
    posn2[1]+w, posn2[0]-w, posn2[0]+w)
wz = zoomSz/2
subF1 = rawFlat[posn1[1]-wz:posn1[1]+wz, posn1[0]-wz:posn1[0]+wz]
subF2 = rawFlat[posn2[1]-wz:posn2[1]+wz, posn2[0]-wz:posn2[0]+wz]


plt.figure(1)
plt.clf()
plt.imshow(subflat1)
plt.figure(2)
plt.clf()
plt.imshow(subflat2)
plt.pause(0.001)

# # Again checking noise
# wz = zoomSz/2
# # subF1 = f1[posn1[1]-wz:posn1[1]+wz, posn1[0]-wz:posn1[0]+wz]
# # subF2 = f2[posn1[1]-wz:posn1[1]+wz, posn1[0]-wz:posn1[0]+wz]
# diff=subF1-subF2
# print np.std(diff)

#np.savez('flats750Aug2015.npz',flatC1=subflat1, flatC2=subflat2 )


if makeManually:
    # Just manually make Make masks
    # Set it as a region between two diagonal lines
    inds = np.arange(0, imsz)
    mask1 = np.ones((imsz,imsz))
    # plt.clf()
    # plt.imshow(subflat1)
    m1_1 = 1.5
    b1_1 = 110.
    # y = m1_1*inds + b1_1
    # plt.plot(y, inds)
    m1_2 = 1.1
    b1_2 = -180
    # y = m1_2*inds + b1_2
    # plt.plot(y, inds)

    yv, xv = np.meshgrid(inds, inds)
    mask1[np.where(yv > m1_1*xv + b1_1)] = 0
    mask1[np.where(yv < m1_2*xv + b1_2)] = 0
    subflat1_mskd  = copy(subflat1)
    subflat1_mskd = subflat1_mskd * mask1

    plt.clf()
    plt.imshow(subflat1_mskd)

    mask2 = np.ones((imsz,imsz))
    plt.clf()
    plt.imshow(subflat2)
    m2_1 = 1.1
    b2_1 = 155.
    # y = m2_1*inds + b2_1
    # plt.plot(y, inds)
    m2_2 = 0.6
    b2_2 = -70
    # y = m2_2*inds + b2_2
    # plt.plot(y, inds)

    mask2[np.where(yv > m2_1*xv + b2_1)] = 0
    mask2[np.where(yv < m2_2*xv + b2_2)] = 0
    subflat2_mskd  = copy(subflat2)
    subflat2_mskd = subflat2_mskd * mask2

    plt.clf()
    plt.imshow(subflat2_mskd)

    np.savez('manualMasksNoc2016.npz',mask1=mask1, mask2=mask2 )