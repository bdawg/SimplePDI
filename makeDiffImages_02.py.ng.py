
import numpy as np
from scipy import ndimage
from astropy.io import fits
import matplotlib.pyplot as plt
import polzimtools as plz



#bsRot = 74.74 #Differential rotation of beamsplitter, in degrees
bsRot = 89.93 # Why is it not 74.74, like in make_tmpl????????

doMask = False
maskRad = 10 # Radius of central mask (px)
maskVal = 1  # Value to set central masked region to

clim=[-5, 5]
clim=[-100, 100]

# Make the various difference and ratio images
loadFilename = 'allSummedImsCube_176files__cube_HD169142_01_20160721_675-50_EmptySlot_.npy'
loadFilename = 'allSummedImsCube2_cube_HD169142_01_20160721_675-50_EmptySlot_CentroidLLim500.npy'
loadFilename = 'allSummedImsCube_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
#loadFilename = 'allSummedImsCube_cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_.npy'
#loadFilename = 'allSummedImsCube_cube_HD19820_01_20160918_750-50_EmptySlot_.npz'

deRotate = True

# Only use a subset of images?
imRange = [20, 164]

# Reload cube if allSummedIms already exists? (=False useful if in interactive mode)
#loadImages = True



###################################


def makeDiffIms(filename, doMask = False, maskRad = 2., maskVal = 1., clim = [], deRotate = True,
                imRange = [], reload = False):

    npzfile = np.load(loadFilename)
    allSummedImsOrig = npzfile['arr_0']
    allPAsIn = npzfile['arr_1']
    del (npzfile)
    imSz = allSummedImsOrig.shape[0]
    nIms = allSummedImsOrig.shape[2]
    if len(imRange) > 0:
        imInd = range(imRange[0], imRange[1])
    else:
        imInd = range(nIms)
    paInds = range(0,len(allPAs),4)
    allPAs = allPAsIn[paInds]


    if deRotate:
        allSummedIms = np.zeros_like(allSummedImsOrig)
        # Approximate set of 4 HWP files as having single PA (the average over files)
        pasPerInd = []
        curPaSet = []
        for ind in imInd:
            print "ind: %d" % ind
            for h in range(4):
                #pa = allPAs[count]
                curPaSet.append(allPAs(count))

                for c in range(2):
                    for l in range(2):
                        curImIn = allSummedImsOrig[:, :, ind, h, c, l]
                        # curImIn = allSummedImsOrig[:, :, ind, h, c, l] - h0c0l0med
                        if c == 0:
                            # Rotate channel0 images
                            curImIn = ndimage.interpolation.rotate(np.fliplr(curImIn), bsRot,
                                                                   mode='nearest', reshape=False)
                        allSummedIms[:, :, ind, h, c, l] = \
                            ndimage.interpolation.rotate(curImIn, -pa, mode='nearest', reshape=False)


            pasPerInd.append(pa)
        pasPerInd = np.asarray(pasPerInd)
    else:












