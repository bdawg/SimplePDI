

fileExtn = '.fits'

dataPath = '/Volumes/silo4/snert/VAMPIRESData_201607/Analysis/HD169142_21072016/'
filePref = 'cube_HD169142_01_20160721_675-50_EmptySlot_'

#dataPath = '/Volumes/silo4/snert/VAMPIRESData201607/Analysis/HD317501--cal--_21072016/'
#filePref = 'cube_HD317501_01_20160721_675-50_EmptySlot_'

dataPath = '/Volumes/silo4/snert/VAMPIRESData_201609/Analysis/ABAur/'
filePref = 'cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_'

dataPath = '/Volumes/silo4/snert/VAMPIRESData_201609/Analysis/MWC758/'
filePref = 'cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_'


cubeInfoFile = 'cubeinfoOct2016.idlvar'

saveFilename = 'allSummedImsCube_'+filePref

# Assumes HWP angles go as 0, 22.5, 45, 67.5 
startFileNum = 0 
nFiles = 164 *4 # NB actual files, not filenums. So 4x num raw files
nFiles = 172 *4 # NB actual files, not filenums. So 4x num raw files

showImage = False

centroidThresh = 5000 #Values below this clipped to zero for centroiding
#centroidThresh = 2000
###########################

import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt
from copy import copy


class readCubeInfo:
    def __init__(self, cubeInfoFileString):
        # Get useful metadata from cubeinfo file
        cubeinfoObj = io.readsav(cubeInfoFileString, python_dict=False, verbose=False)
        self.UTCs = cubeinfoObj.olog.utc[0]
        self.filters = cubeinfoObj.olog.filter[0]
        self.ras = cubeinfoObj.olog.ra[0]
        self.decs = cubeinfoObj.olog.dec[0]
        self.mask = cubeinfoObj.olog.mask[0]
        self.adate = cubeinfoObj.olog.adate[0]
        self.emgains = cubeinfoObj.olog.emgain[0]
        #self.mffile = cubeinfoObj.plog.mf_file[0]
        self.pkflux = cubeinfoObj.framestats.pkflx[0]
        self.totflux = cubeinfoObj.framestats.totflx[0]
        self.cubename = cubeinfoObj.olog.cube_fname[0][0]
        self.pa = cubeinfoObj.olog.pa[0]
        del (cubeinfoObj)



# Read a FITS file to get sizes
curFilenumStr = '%05d' % startFileNum
curBSStr = '_1'
curLCVRStr = '_A'
curFilename = dataPath + filePref + curFilenumStr + curBSStr + curLCVRStr \
    + fileExtn
hdulist = fits.open(curFilename)
curHDU = hdulist[0]
curCube = np.transpose(curHDU.data)
nFrms = curCube.shape[2]
dim = curCube.shape[0]

# Get PAs
cubeInfoFileString = dataPath + cubeInfoFile
cubeInfo = readCubeInfo(cubeInfoFileString)
pas = cubeInfo.pa


# Indexes are [:, :, file, HWP, Channel, LCVRstate]
allSummedIms = np.zeros([dim,dim,nFiles/4, 4, 2, 2])
curFileNum = startFileNum
hwpState = 0 # This counts through 4 positions
curHWPSet = 0 # Increments for each new set of HWP posns

for f in range(0, nFiles/4):
    curFileNum = f
    
    for c in range(0,2):
        curBSChan = c 
          
        for l in range(0,2):
            curLCVRState = l

            # Generate current filename
            if curBSChan == 0:
                curBSStr = '_1'
            else:
                curBSStr = '_2'
                
            if curLCVRState == 0:
                curLCVRStr = '_A'
            else:
                curLCVRStr = '_B' 
                
            curFilenumStr = '%05d' % curFileNum              
            curFilename = dataPath + filePref + curFilenumStr + curBSStr + curLCVRStr \
                + fileExtn
            print 'Reading file %s' % curFilename
            print 'HWPSet %d -- HWP State %d, Channel %d, LCVR %d' % \
                (curHWPSet, hwpState, curBSChan, curLCVRState)
            hdulist = fits.open(curFilename)
            curHDU = hdulist[0]
            curCube = np.transpose(curHDU.data)
    
            # Show summed image
            curCubeSummed = np.mean(curCube, axis=2)
            if showImage:
                fig1 = plt.figure(1)
                plt.imshow(curCubeSummed, interpolation='nearest')
                plt.pause(0.1)
            print 'Max val in summed UNshifted image: %f' % np.max(curCubeSummed.ravel())

            # Make cube of shifted images
            curCubeShifted = np.zeros([dim, dim, nFrms])
            # Centre the image
            for i in range(0, nFrms):

                llim = centroidThresh
                curIm = copy(curCube[:,:,i])
                #print np.min(curIm.min())
                curIm[np.where(curIm <= llim)] = 0
                curCent = ndimage.center_of_mass(curIm)
                #curCent = np.where(curCube[:,:,i] == curCube[:,:,i].max())
                #print curCent

                if showImage:
                    plt.clf()
                    plt.subplot(121)
                    plt.imshow(curCube[:,:,i], interpolation='nearest')
                    plt.plot(curCent[1], curCent[0], 'rx', ms=20, mew=1)
                    plt.subplot(122)
                    plt.imshow(curIm, interpolation='nearest')
                    plt.plot(curCent[1], curCent[0], 'rx', ms=20, mew=1)
                    plt.pause(0.001)

                imCent = dim/2
                rowOffset = imCent-curCent[0]
                colOffset = imCent-curCent[1]

                curCubeShifted[:,:,i] = ndimage.interpolation.shift(curCube[:,:,i],
                    [rowOffset, colOffset])
                               
            curCubeShiftedSummed = np.mean(curCubeShifted, axis=2)           
            curHDU = fits.PrimaryHDU()
            curHDU.data = curCubeShiftedSummed
            curHDU.writeto('autosave.fits', clobber=True)

            if showImage:
                plt.clf()
                plt.imshow(curCubeShiftedSummed, interpolation='nearest')
                plt.pause(0.1)
            print 'Max val in summed shifted image:   %f' % np.max(curCubeShiftedSummed.ravel())
            print ' '
            allSummedIms[:,:,f,hwpState,c,l] = curCubeShiftedSummed
            
    # Increment HWP state
    hwpState = hwpState + 1
    if hwpState == 4:
        hwpState = 0
        curHWPSet = curHWPSet + 1
    
    
# Save the big cube
#np.save(saveFilename,allSummedIms)
np.savez(saveFilename, allSummedIms, pas)



# Make individual summed images
# Indexes are [:, :, file, HWP, Channel, LCVRstate]
h0c0l0 = np.mean(allSummedIms[:, :, :, 0, 0, 0], axis=2)
h0c0l1 = np.mean(allSummedIms[:, :, :, 0, 0, 1], axis=2)
h0c1l0 = np.mean(allSummedIms[:, :, :, 0, 1, 0], axis=2)
h0c1l1 = np.mean(allSummedIms[:, :, :, 0, 1, 1], axis=2)
h1c0l0 = np.mean(allSummedIms[:, :, :, 1, 0, 0], axis=2)
h1c0l1 = np.mean(allSummedIms[:, :, :, 1, 0, 1], axis=2)
h1c1l0 = np.mean(allSummedIms[:, :, :, 1, 1, 0], axis=2)
h1c1l1 = np.mean(allSummedIms[:, :, :, 1, 1, 1], axis=2)
h2c0l0 = np.mean(allSummedIms[:, :, :, 2, 0, 0], axis=2)
h2c0l1 = np.mean(allSummedIms[:, :, :, 2, 0, 1], axis=2)
h2c1l0 = np.mean(allSummedIms[:, :, :, 2, 1, 0], axis=2)
h2c1l1 = np.mean(allSummedIms[:, :, :, 2, 1, 1], axis=2)
h3c0l0 = np.mean(allSummedIms[:, :, :, 3, 0, 0], axis=2)
h3c0l1 = np.mean(allSummedIms[:, :, :, 3, 0, 1], axis=2)
h3c1l0 = np.mean(allSummedIms[:, :, :, 3, 1, 0], axis=2)
h3c1l1 = np.mean(allSummedIms[:, :, :, 3, 1, 1], axis=2)

for h in range(0,4):
    for c in range(0,2):
        for l in range(0,2):
            idString = 'HWP%d_C%d_L%d' % (h, c, l) 
            fitsFilename = 'summedState_' + filePref + idString
            print 'Saving FITS file ' + fitsFilename
            curIm = np.mean(allSummedIms[:, :, :, h, c, l], axis=2)
            curHDU = fits.PrimaryHDU()
            curHDU.data = curIm
            curHDU.writeto(fitsFilename)











