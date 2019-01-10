# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 18:45:49 2016

@author: bnorris
"""
import numpy as np
from scipy import ndimage
from astropy.io import fits
import matplotlib.pyplot as plt

#bsRot = 74.74 #Differential rotation of beamsplitter, in degrees
bsRot = 89.93 #Why is it not 74.74, like in make_tmpl????????


doMask = False
maskRad = 10 #Radius of central mask (px)
maskVal = 1  #Value to set central masked region to 

clim=[-5, 5]
clim=[-100, 100]

# Make the various difference and ratio images
loadFilename = 'allSummedImsCube_176files__cube_HD169142_01_20160721_675-50_EmptySlot_.npy'
loadFilename = 'allSummedImsCube2_cube_HD169142_01_20160721_675-50_EmptySlot_CentroidLLim500.npy'
loadFilename = 'allSummedImsCube_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
#loadFilename = 'allSummedImsCube_cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_.npy'

doDerotate = True

# Only use a subset of images?
imRange = [20, 164]

# Reload cube if allSummedIms already exists?
#reload = False
###################################

# try: allSummedIms
# except:
#     allSummedIms = np.load(loadFilename)
# else:
#     if reload:
#         allSummedIms = np.load(loadFilename)

npzfile = np.load(loadFilename)
allSummedImsOrig = npzfile['arr_0']
allPAs = npzfile['arr_1']
del(npzfile)

imSz = allSummedImsOrig.shape[0]
nIms = allSummedImsOrig.shape[2]

try: imRange
except NameError: imInd = range(nIms)
else: imInd = range(imRange[0], imRange[1])



#h0c0l0med = np.mean(allSummedImsOrig[:, :, imInd, 0, 0, 0], axis=2)


if doDerotate:
    # Derotate images
    allSummedIms = np.zeros_like(allSummedImsOrig)
    count = 0
    for ind in imInd:
        print "ind: %d" % ind
        for h in range(4):
            pa = allPAs[count]
            for c in range(2):
                for l in range(2):
                    curImIn = allSummedImsOrig[:,:,ind,h,c,l]
                    #curImIn = allSummedImsOrig[:, :, ind, h, c, l] - h0c0l0med
                    if c == 0:
                        # Rotate channel0 images
                        curImIn = ndimage.interpolation.rotate(np.fliplr(curImIn), bsRot, mode='nearest', reshape=False)
                    allSummedIms[:,:,ind,h,c,l] = \
                        ndimage.interpolation.rotate(curImIn,-pa, mode='nearest', reshape=False)

            count += 1
else:
    allSummedIms = allSummedImsOrig



# Make individual summed images
# Indexes are [:, :, file, HWP, Channel, LCVRstate]
h0c0l0 = np.mean(allSummedIms[:, :, imInd, 0, 0, 0], axis=2)
h0c0l1 = np.mean(allSummedIms[:, :, imInd, 0, 0, 1], axis=2)
h0c1l0 = np.mean(allSummedIms[:, :, imInd, 0, 1, 0], axis=2)
h0c1l1 = np.mean(allSummedIms[:, :, imInd, 0, 1, 1], axis=2)
h1c0l0 = np.mean(allSummedIms[:, :, imInd, 1, 0, 0], axis=2)
h1c0l1 = np.mean(allSummedIms[:, :, imInd, 1, 0, 1], axis=2)
h1c1l0 = np.mean(allSummedIms[:, :, imInd, 1, 1, 0], axis=2)
h1c1l1 = np.mean(allSummedIms[:, :, imInd, 1, 1, 1], axis=2)
h2c0l0 = np.mean(allSummedIms[:, :, imInd, 2, 0, 0], axis=2)
h2c0l1 = np.mean(allSummedIms[:, :, imInd, 2, 0, 1], axis=2)
h2c1l0 = np.mean(allSummedIms[:, :, imInd, 2, 1, 0], axis=2)
h2c1l1 = np.mean(allSummedIms[:, :, imInd, 2, 1, 1], axis=2)
h3c0l0 = np.mean(allSummedIms[:, :, imInd, 3, 0, 0], axis=2)
h3c0l1 = np.mean(allSummedIms[:, :, imInd, 3, 0, 1], axis=2)
h3c1l0 = np.mean(allSummedIms[:, :, imInd, 3, 1, 0], axis=2)
h3c1l1 = np.mean(allSummedIms[:, :, imInd, 3, 1, 1], axis=2)

# # Rotate channel0 images
# h0c0l0 = ndimage.interpolation.rotate(np.fliplr(h0c0l0), bsRot, reshape=False)
# h0c0l1 = ndimage.interpolation.rotate(np.fliplr(h0c0l1), bsRot, reshape=False)
# h1c0l0 = ndimage.interpolation.rotate(np.fliplr(h1c0l0), bsRot, reshape=False)
# h1c0l1 = ndimage.interpolation.rotate(np.fliplr(h1c0l1), bsRot, reshape=False)
# h2c0l0 = ndimage.interpolation.rotate(np.fliplr(h2c0l0), bsRot, reshape=False)
# h2c0l1 = ndimage.interpolation.rotate(np.fliplr(h2c0l1), bsRot, reshape=False)
# h3c0l0 = ndimage.interpolation.rotate(np.fliplr(h3c0l0), bsRot, reshape=False)
# h3c0l1 = ndimage.interpolation.rotate(np.fliplr(h3c0l1), bsRot, reshape=False)

# Apply mask
if doMask:
    maskPosn = imSz/2
    yGrid, xGrid = np.ogrid[-maskPosn:imSz-maskPosn, -maskPosn:imSz-maskPosn]
    maskInds = xGrid**2 + yGrid**2 <= maskRad**2
    h0c0l0[maskInds] = maskVal
    h0c0l1[maskInds] = maskVal
    h0c1l0[maskInds] = maskVal
    h0c1l1[maskInds] = maskVal
    h1c0l0[maskInds] = maskVal
    h1c0l1[maskInds] = maskVal
    h1c1l0[maskInds] = maskVal
    h1c1l1[maskInds] = maskVal
    h2c0l0[maskInds] = maskVal
    h2c0l1[maskInds] = maskVal
    h2c1l0[maskInds] = maskVal
    h2c1l1[maskInds] = maskVal
    h3c0l0[maskInds] = maskVal
    h3c0l1[maskInds] = maskVal
    h3c1l0[maskInds] = maskVal
    h3c1l1[maskInds] = maskVal



allFlattenedIms = np.asarray([h0c0l0, h0c0l1, h0c1l0, h0c1l1, h1c0l0, h1c0l1, 
                              h1c1l0, h1c1l1, h2c0l0, h2c0l1, h2c1l0, h2c1l1,
                              h3c0l0, h3c0l1, h3c1l0, h3c1l1])
                              
idStrings = ['h0c0l0', 'h0c0l1', 'h0c1l0', 'h0c1l1', 'h1c0l0', 'h1c0l1', 
             'h1c1l0', 'h1c1l1', 'h2c0l0', 'h2c0l1', 'h2c1l0', 'h2c1l1',
             'h3c0l0', 'h3c0l1', 'h3c1l0', 'h3c1l1']                              
                              
plt.clf()
for k in range(16):
    im = allFlattenedIms[k,:,:]
    plt.subplot(4,4,k+1)
    plt.imshow(im, clim=clim, interpolation='Nearest')
    plt.title(idStrings[k])


#h0c0l0 = np.clip(h0c0l0,-5,50)
#h0c0l1 = np.clip(h0c0l1,-5,50)
#h0c1l0 = np.clip(h0c1l0,-5,50)
#h0c1l1 = np.clip(h0c1l1,-5,50)
#
#im1 = h0c0l0-h0c0l1
#im2 = h0c1l0-h0c1l1
#dd = im2-im1


"""
Here, use the double-ratio method
e.g. Quanz, et al. 2013 (HD 169142)
"""

# Do with no HWP first
#RQ1 = np.sqrt( (h0c0l0 / h0c1l0) / (h0c0l1 / h0c1l1) )
#pQ1 = (RQ1-1)/(RQ1+1)
#I1 = h0c0l0 + h0c1l1
#Q1 = pQ1 / I1
#
#RQ1 = np.sqrt( (h0c0l0 / h0c1l0) / (h2c0l0 / h2c1l0) )
#pQ1 = (RQ1-1)/(RQ1+1)
#I1 = h0c0l0 + h0c1l1
#Q1 = pQ1 / I1
#
#RQ1 = np.sqrt( (h1c0l0 / h1c1l0) / (h3c0l0 / h3c1l0) )
#pQ1 = (RQ1-1)/(RQ1+1)
#I1 = h0c0l0 + h0c1l1
#Q1 = pQ1 / I1
#
#RQ1 = np.sqrt( (h0c1l0 / h0c1l1) / (h2c1l0 / h2c1l1) )
#pQ1 = (RQ1-1)/(RQ1+1)
#I1 = h0c0l0 + h0c1l1
#Q1 = pQ1 / I1
#
#im=Q1
#im[np.isnan(im)] = 1



plusQ1 = h0c0l0 - h0c1l0
minusQ1 = h0c0l1 - h0c1l1
Q1 = (plusQ1 - minusQ1)/2
#I0 = (h0c0l0 + h0c1l1)/2
I0 = (h0c0l0 + h0c1l0 + h0c0l1 + h0c1l1)/4

plusQ2 = h2c0l0 - h2c1l0
minusQ2 = h2c0l1 - h2c1l1
Q2 = (plusQ2 - minusQ2)/2
I2 = (h2c0l0 + h2c1l0 + h2c0l1 + h2c1l1)/4

Q = Q1 - Q2
I_02 = (I0 + I2)/2

pQ=Q/I_02


plusU1 = h1c0l0 - h1c1l0
minusU1 = h1c0l1 - h1c1l1
U1 = (plusU1 - minusU1)/2
#I0 = (h0c0l0 + h0c1l1)/2
I1 = (h1c0l0 + h1c1l0 + h1c0l1 + h1c1l1)/4

plusU2 = h3c0l0 - h3c1l0
minusU2 = h3c0l1 - h3c1l1
U2 = (plusU2 - minusU2)/2
I3 = (h3c0l0 + h3c1l0 + h3c0l1 + h3c1l1)/4

U = U1 - U2
I_13 = (I1 + I3)/2

pU=U/I_13

p = np.sqrt(pQ**2 + pU**2)



im=Q
#im = np.clip(im,-0.1,0.1)

#im=h0c0l0
#im=h0c0l0

im = pQ

curHDU = fits.PrimaryHDU()
curHDU.data = im
curHDU.writeto('out.fits', clobber=True)


print 'Done.'