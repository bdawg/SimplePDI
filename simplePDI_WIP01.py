# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 16:46:37 2016

@author: bnorris
"""

import numpy as np
from scipy import ndimage
from astropy.io import fits
import matplotlib.pyplot as plt
import time

dataPath = '/Volumes/silo4/snert/VAMPIRESData201607/Analysis/HD169142_21072016/'
filePref = 'cube_HD169142_01_20160721_675-50_EmptySlot_'
fileExtn = '.fits'

#dataPath = '/Volumes/silo4/snert/VAMPIRESData201607/Analysis/HD317501--cal--_21072016/'
#filePref = 'cube_HD317501_01_20160721_675-50_EmptySlot_'

curFileNum = 68
curBSChan = 1 # 1 or 2
curLCVRState = 1 # 1 or 2


# Generate current filename
if curBSChan == 1:
    curBSStr = '_1'
else:
    curBSStr = '_2'
    
if curLCVRState == 1:
    curLCVRStr = '_A'
else:
    curLCVRStr = '_B' 
    
curFilenumStr = '%05d' % curFileNum
    
curFilename = dataPath + filePref + curFilenumStr + curBSStr + curLCVRStr \
    + fileExtn

#curFilename = '/Volumes/silo4/snert/VAMPIRESData201607/20160721/HD169142_01_20160721_675-50_EmptySlot_04.fits'


# Read FITS file
hdulist = fits.open(curFilename)
curHDU = hdulist[0]
curCube = np.transpose(curHDU.data)
nFrms = curCube.shape[2]
dim = curCube.shape[0]

curCubeSummed = np.sum(curCube, axis=2)

# Show summed image
fig1 = plt.figure(1)
plt.imshow(curCubeSummed, interpolation='nearest')
plt.pause(0.1)
print 'Max val in summed UNshifted image: %f' % np.max(curCubeSummed.ravel())


curCubeShifted = np.zeros([dim, dim, nFrms])
# Centre the image
for i in range(0, nFrms):
    #curCent = ndimage.center_of_mass(curCube[:,:,i])   
    curCent = np.where(curCube[:,:,i] == curCube[:,:,i].max())
 
#    print curCent
#    plt.clf()
#    plt.imshow(curCube[:,:,i], interpolation='nearest')
#    plt.plot(curCent[1], curCent[0], 'ro')
#    plt.pause(0.001)
    
    imCent = dim/2
    rowOffset = imCent-curCent[0]
    colOffset = imCent-curCent[1]
    
    curCubeShifted[:,:,i] = ndimage.interpolation.shift(curCube[:,:,i], 
        [rowOffset, colOffset])
        
    print i
    
    
curCubeShiftedSummed = np.sum(curCubeShifted, axis=2)

curHDU = fits.PrimaryHDU()
curHDU.data = curCubeShiftedSummed
curHDU.writeto('autosave.fits', clobber=True)


plt.clf()
plt.imshow(curCubeShiftedSummed, interpolation='nearest')
print 'Max val in summed shifted image: %f' % np.max(curCubeShiftedSummed.ravel())
