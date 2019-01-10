"""
TODO:
- Try changing order of divisions (and excluding some). Use some matrix approach??
- Try flat-fielding. (And dark?)
- Try more averaging before dividing
- Try rotating more finely in time
"""


import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt
from copy import copy
import polzimtools as plz
import time
infile = 'allSummedImsSaved_abaur01.npz'

#np.seterr(all='raise')

class sortPolzState:
    def __init__(self, dataIn, ind):
        self.h0c0l0 = dataIn[:, :, ind, 0, 0, 0]
        self.h0c0l1 = dataIn[:, :, ind, 0, 0, 1]
        self.h0c1l0 = dataIn[:, :, ind, 0, 1, 0]
        self.h0c1l1 = dataIn[:, :, ind, 0, 1, 1]
        self.h1c0l0 = dataIn[:, :, ind, 1, 0, 0]
        self.h1c0l1 = dataIn[:, :, ind, 1, 0, 1]
        self.h1c1l0 = dataIn[:, :, ind, 1, 1, 0]
        self.h1c1l1 = dataIn[:, :, ind, 1, 1, 1]
        self.h2c0l0 = dataIn[:, :, ind, 2, 0, 0]
        self.h2c0l1 = dataIn[:, :, ind, 2, 0, 1]
        self.h2c1l0 = dataIn[:, :, ind, 2, 1, 0]
        self.h2c1l1 = dataIn[:, :, ind, 2, 1, 1]
        self.h3c0l0 = dataIn[:, :, ind, 3, 0, 0]
        self.h3c0l1 = dataIn[:, :, ind, 3, 0, 1]
        self.h3c1l0 = dataIn[:, :, ind, 3, 1, 0]
        self.h3c1l1 = dataIn[:, :, ind, 3, 1, 1]
        self.combined = np.asarray(
            [self.h0c0l0, self.h0c0l1, self.h0c1l0, self.h0c1l1, self.h1c0l0,
             self.h1c0l1, self.h1c1l0, self.h1c1l1, self.h2c0l0, self.h2c0l1,
             self.h2c1l0, self.h2c1l1, self.h3c0l0, self.h3c0l1, self.h3c1l0, self.h3c1l1])
        self.idStrings = ['h0c0l0', 'h0c0l1', 'h0c1l0', 'h0c1l1', 'h1c0l0', 'h1c0l1',
                          'h1c1l0', 'h1c1l1', 'h2c0l0', 'h2c0l1', 'h2c1l0', 'h2c1l1',
                          'h3c0l0', 'h3c0l1', 'h3c1l0', 'h3c1l1']

class sortPolzStateTesting:
    # This test class changes round the states used (so they're actually mislabelled)
    def __init__(self, dataIn, ind):
        # # This version - swap chan and HWP
        # self.h0c0l0 = dataIn[:, :, ind, 0, 0, 0]
        # self.h0c0l1 = dataIn[:, :, ind, 0, 0, 1]
        # self.h0c1l0 = dataIn[:, :, ind, 2, 0, 0]
        # self.h0c1l1 = dataIn[:, :, ind, 2, 0, 1]
        # self.h2c0l0 = dataIn[:, :, ind, 0, 1, 0]
        # self.h2c0l1 = dataIn[:, :, ind, 0, 1, 1]
        # self.h2c1l0 = dataIn[:, :, ind, 2, 1, 0]
        # self.h2c1l1 = dataIn[:, :, ind, 2, 1, 1]
        #
        # self.h1c0l0 = dataIn[:, :, ind, 1, 0, 0]
        # self.h1c0l1 = dataIn[:, :, ind, 1, 0, 1]
        # self.h1c1l0 = dataIn[:, :, ind, 3, 0, 0]
        # self.h1c1l1 = dataIn[:, :, ind, 3, 0, 1]
        # self.h3c0l0 = dataIn[:, :, ind, 1, 1, 0]
        # self.h3c0l1 = dataIn[:, :, ind, 1, 1, 1]
        # self.h3c1l0 = dataIn[:, :, ind, 3, 1, 0]
        # self.h3c1l1 = dataIn[:, :, ind, 3, 1, 1]

        # This version - swap lcvr and HWP
        self.h0c0l0 = dataIn[:, :, ind, 0, 0, 0]
        self.h0c0l1 = dataIn[:, :, ind, 2, 0, 0]
        self.h0c1l0 = dataIn[:, :, ind, 0, 1, 0]
        self.h0c1l1 = dataIn[:, :, ind, 2, 1, 0]
        self.h2c0l0 = dataIn[:, :, ind, 0, 0, 1]
        self.h2c0l1 = dataIn[:, :, ind, 2, 0, 1]
        self.h2c1l0 = dataIn[:, :, ind, 0, 1, 1]
        self.h2c1l1 = dataIn[:, :, ind, 2, 1, 1]

        self.h1c0l0 = dataIn[:, :, ind, 1, 0, 0]
        self.h1c0l1 = dataIn[:, :, ind, 3, 0, 0]
        self.h1c1l0 = dataIn[:, :, ind, 1, 1, 0]
        self.h1c1l1 = dataIn[:, :, ind, 3, 1, 0]
        self.h3c0l0 = dataIn[:, :, ind, 1, 0, 1]
        self.h3c0l1 = dataIn[:, :, ind, 3, 0, 1]
        self.h3c1l0 = dataIn[:, :, ind, 1, 1, 1]
        self.h3c1l1 = dataIn[:, :, ind, 3, 1, 1]

        self.combined = np.asarray(
            [self.h0c0l0, self.h0c0l1, self.h0c1l0, self.h0c1l1, self.h1c0l0,
             self.h1c0l1, self.h1c1l0, self.h1c1l1, self.h2c0l0, self.h2c0l1,
             self.h2c1l0, self.h2c1l1, self.h3c0l0, self.h3c0l1, self.h3c1l0, self.h3c1l1])
        self.idStrings = ['h0c0l0', 'h0c0l1', 'h0c1l0', 'h0c1l1', 'h1c0l0', 'h1c0l1',
                          'h1c1l0', 'h1c1l1', 'h2c0l0', 'h2c0l1', 'h2c1l0', 'h2c1l1',
                          'h3c0l0', 'h3c0l1', 'h3c1l0', 'h3c1l1']



npzobj = np.load(infile)
allSummedIms = npzobj['allSummedIms']
pasPerInd = npzobj['pasPerInd']
nUsedSets = allSummedIms.shape[2]
#cur = sortPolzState(allSummedIms, 0)
imSz = allSummedIms.shape[0]

allpQIms = []
allpUIms = []

for count in range(nUsedSets):
    cur = sortPolzState(allSummedIms, count)

    try:
        RQ1 = np.sqrt( (cur.h0c0l0 / cur.h0c1l0) / (cur.h0c0l1 / cur.h0c1l1) )
        RQ2 = np.sqrt( (cur.h2c0l0 / cur.h2c1l0) / (cur.h2c0l1 / cur.h2c1l1) )
        RQ = np.sqrt( RQ1 / RQ2 )
        #### RQ = RQ2 #### ignore HWP
        pQ = (RQ - 1) / (RQ + 1)

        RU1 = np.sqrt( (cur.h1c0l0 / cur.h1c1l0) / (cur.h1c0l1 / cur.h1c1l1) )
        RU2 = np.sqrt( (cur.h3c0l0 / cur.h3c1l0) / (cur.h3c0l1 / cur.h3c1l1) )
        RU = np.sqrt( RU1 / RU2 )
        #### RU = RU2 #### ignore HWP
        pU = (RU - 1) / (RU + 1)
        #cur_p_frac = np.sqrt(pQ**2 + pU**2)

        # Make dummy Stokes vector so can reuse the optimised rotImPolz
        curIm = np.zeros((imSz, imSz, 4))
        curIm[:, :, 1] = pQ
        curIm[:, :, 2] = pU

        curIm = plz.rotImPolz(curIm, pasPerInd[count])
    except:
        print "Error caught!"
    else:
        allpQIms.append(curIm[:, :, 1])
        allpUIms.append(curIm[:, :, 2])
    print count

#summedpQIms = np.asarray
summedpQ = np.asarray(allpQIms).mean(axis=0)
summedpU = np.asarray(allpUIms).mean(axis=0)
summed_p = np.sqrt(summedpQ**2 + summedpU**2)

#fits.writeto('testorder_normal.fits',summed_p)
plt.imshow(summed_p)
plt.pause(0.001)


# for count in range(nUsedSets):
#     cur = sortPolzStateTesting(allSummedIms, count)
#
#     RQ1 = np.sqrt( (cur.h0c0l0 / cur.h0c1l0) / (cur.h0c0l1 / cur.h0c1l1) )
#     RQ2 = np.sqrt( (cur.h2c0l0 / cur.h2c1l0) / (cur.h2c0l1 / cur.h2c1l1) )
#     RQ = np.sqrt( RQ1 / RQ2 )
#     #RQ = RQ2  #### ignore "HWP"
#     pQ = (RQ - 1) / (RQ + 1)
#
#     RU1 = np.sqrt( (cur.h1c0l0 / cur.h1c1l0) / (cur.h1c0l1 / cur.h1c1l1) )
#     RU2 = np.sqrt( (cur.h3c0l0 / cur.h3c1l0) / (cur.h3c0l1 / cur.h3c1l1) )
#     RU = np.sqrt( RU1 / RU2 )
#     #RU = RU2  #### ignore "HWP"
#     pU = (RU - 1) / (RU + 1)
#     #cur_p_frac = np.sqrt(pQ**2 + pU**2)
#
#     # Make dummy Stokes vector so can reuse the optimised rotImPolz
#     curIm = np.zeros((imSz, imSz, 4))
#     curIm[:, :, 1] = pQ
#     curIm[:, :, 2] = pU
#
#     curIm = plz.rotImPolz(curIm, pasPerInd[count])
#
#     allpQIms.append(curIm[:, :, 1])
#     allpUIms.append(curIm[:, :, 2])
#
#     print count
#
# #summedpQIms = np.asarray
# summedpQ = np.asarray(allpQIms).mean(axis=0)
# summedpU = np.asarray(allpUIms).mean(axis=0)
# summed_pTesting = np.sqrt(summedpQ**2 + summedpU**2)
#
# #fits.writeto('testorder.fits',summed_pTesting)
# plt.imshow(summed_pTesting)
# plt.pause(0.001)












