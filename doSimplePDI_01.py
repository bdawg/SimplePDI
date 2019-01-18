"""
Script to run simplePDI with various tasks
"""

from simplePDI import *


#########
# Generic options - may be overriden by defaults
showAllImages = False
showPlots = True
saveAllSummedIms = False

centroidThresh = 5000 #Values below this clipped to zero for centroiding
#centroidThresh = 2000

pixScale = 6.5 # Just for scale bar



######################################################################################################
# Usage for reading *raw* files (no IDL pre-cubing)
IDLPreCubed = False






######################################################################################################
# Settings for reading the IDL FITS files (cubed from IDL)
fileExtn = '.fits'

IDLPreCubed = True # True for IDL pipeline cubes, False for raw camera files (only for two-cam data)

# dataPath = '/import/silo4/snert/VAMPIRESData_201607/Analysis/HD169142_21072016/'
# filePref = 'cube_HD169142_01_20160721_675-50_EmptySlot_'

#dataPath = '/Volumes/silo4/snert/VAMPIRESData201607/Analysis/HD317501--cal--_21072016/'
#filePref = 'cube_HD317501_01_20160721_675-50_EmptySlot_'

#dataPath = '/import/silo4/snert/VAMPIRESData_201609/Analysis/ABAur/'
dataPath = "/Users/bnorris/DontBackup/simplePDIdata/ABAur/"
filePref = 'cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_'
nSubFiles = 164 *4 # NB actual files, not filenums. So 4x num raw files
cubeInfoFile = 'cubeinfoOct2016.idlvar'

# dataPath = '/import/silo4/snert/VAMPIRESData_201609/Analysis/MWC758/'
# filePref = 'cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_'
# nSubFiles = 172 *4 # NB actual files, not filenums. So 4x num raw files

# dataPath = "/Users/bnorris/DontBackup/simplePDIdata/ABAur_201610/"
# filePref = 'cube_ABAur_01-02Combined_20161015_750-50_EmptySlot_'
# nSubFiles = 116 *4 # NB actual files, not filenums. So 4x num raw files
# cubeInfoFile = 'cubeinfoFeb2017.idlvar'

saveFilePref = 'allSummedImsCube_20170227_maxmeth'
saveFilePref = 'WIP'

# Assumes HWP angles go as 0, 22.5, 45, 67.5
startFileNum = 0
# nSubFiles = 164 *4 # NB actual files, not filenums. So 4x num raw files
# nSubFiles = 172 *4 # NB actual files, not filenums. So 4x num raw files






######################################################################################################
# Settings for making the differential images

twoCamMode = True # If true, just reflect chan2 horizontally, not rotate.

#bsRot = 74.74 #Differential rotation of beamsplitter, in degrees
bsRot = 89.93 #Why is it not 74.74, like in make_tmpl????????


doMask = False
maskRad = 10 #Radius of central mask (px)
maskVal = 1  #Value to set central masked region to

flatFileC1 = ''
flatFileC2 = ''

clim=[-5, 5]
clim=[-100, 100]

dataPath = './'
dataPath = "/Users/bnorris/DontBackup/simplePDIdata/"
dataPath = "/Volumes/RedDisk/VAMPIRES_WorkingData/simplePDIdata/"
#dataPath = '/Volumes/BN_DATA_OSX/VAMPIRES_WorkingData/MWC758/'
# dataPath = '/Users/bnorris/DataAnalysis/VAMPIRES_DataAnalysis/SimplePDI/'
# dataPath = '/Users/bnorris/DontBackup/simplePDIdata/'

# loadFilename = 'allSummedImsCube_176files__cube_HD169142_01_20160721_675-50_EmptySlot_.npy'
# loadFilename = 'allSummedImsCube2_cube_HD169142_01_20160721_675-50_EmptySlot_CentroidLLim500.npy'
# loadFilename = 'allSummedImsCube_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
#loadFilename = 'allSummedImsCube_cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_.npy'
#loadFilename = 'allSummedImsCube_cube_HD19820_01_20160918_750-50_EmptySlot_.npz'
# loadFilename = 'allSummedImsCubeB_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
# loadFilename = 'allSummedImsCube_maxMeth_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
loadFilename = 'allSummedImsCube_centMeth_cube_ABAur_01-04Combined_20160918_750-50_EmptySlot_.npz'
# loadFilename = 'allSummedImsCube_cenMeth_cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_.npz'
# loadFilename = 'allSummedImsCube_maxMeth_cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_.npz'
#loadFilename = dataPath+loadFilename
#loadFilename = 'allSummedImsCube_maxMeth_cube_AUMic20161015_01_Combined_20161015_750-50_EmptySlot_.npz'
loadFilename = 'allSummedImsCube_20170227cube_ABAur_01-02Combined_20161015_750-50_EmptySlot_.npz'
# loadFilename ='allSummedImsCube_20170227_maxmethcube_ABAur_01-02Combined_20161015_750-50_EmptySlot_.npz'
# loadFilename= 'allSummedImsCube_201702ABAur_01-03Combined_20170116_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_201701_332acqs_maxmethABAur_01-03Combined_20170116_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_201701_332acqs_COMmethABAur_01-03Combined_20170116_750-50_EmptySlot_0.npz'
loadFilename = 'allSummedImsCube_ABAur201612_52acqs_COMmethABAur_02_20161215_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_proc201709_COMmeth_ABAur_01_20170313_750-50_EmptySlot_0.npz'

## loadFilename = 'allSummedImsCube_201709_COMmeth_omiCet_01_20170912_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_201710_68acqs_COMmeth_subtDark2_omiCet_01_20170912_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_201710_68acqs_COMmeth_noBGsubt_omiCet_01_20170912_750-50_EmptySlot_0.npz'

# loadFilename = 'allSummedImsCube_201709_COMmeth_PREFLIP_omiCet_01_20170912_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_201711_COMmeth_subtdark_63Cet_01_20170912_750-50_EmptySlot_0.npz'
###loadFilename = 'allSummedImsCube_201704_128acqs_COMmethHD141569_02Combined_20170418_750-50_EmptySlot_0.npz'


dataPath = './'
dataPath = '../SimplePDI_DATA/'
loadFilename = 'allSummedImsCube_MAXmeth_darksubt_ABAur_01_20181017_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_MAXmeth_darksubt_ABAur_02_20181017_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_MAXmeth_darksubt_ABAur_03__RESTARTED___20181017_750-50_EmptySlot_0.npz'
loadFilename = 'allSummedImsCube_MAXmeth_darksubt_ABAur_01_20181018_750-50_EmptySlot_0.npz'

doDerotate = True

# Only use a subset of images? This is in *HWP set* numbering
# imRange = [5, 41]
#imRange = [0, 29]

# Reload cube if allSummedIms already exists?
#reload = False





################################################################################################################
"""
Do stuff:
"""

calMat = np.asarray([[1., 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])

# Includes no 90 deg rotation
calMat = np.asarray([[1.0653, 0.0457, 0.0084, 0],
                     [0.0228, 0.7153, -0.0355, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])
# No 90 deg, but included 3rd row before inv
calMat = np.asarray([[1.0655, 0.0458, 0.0060, 0],
                     [0.0220, 0.7147, -0.0254, 0],
                     [0.0233, 0.0163, 0.7139, 0],
                     [0, 0, 0, 1]])

# # Includes 90 deg rotation
# calMat = np.asarray([[1.0624, -0.0455, 0.0129, 0],
#                      [-0.0228, -0.7133, 0.0354, 0],
#                      [0, 0, 1, 0],
#                      [0, 0, 0, 1]])


#
p = pdiImages()

# p.processRawFiles(dataPath, filePref, cubeInfoFile, nSubFiles, startFileNum=startFileNum,
#                   saveFilePref = saveFilePref, centroidThresh = centroidThresh, showImage=False,
#                   method='max', twoCamMode = False, IDLPreCubed = True)
# p.loadCube(dataPath, loadFilename)
p.loadCube(dataPath, loadFilename, twoCamMode = twoCamMode)

# p.makeDiffIms(imRange=imRange, showPlots=False, twoCamMode=twoCamMode, deRotate=doDerotate)
p.makeDiffIms(showPlots=False, twoCamMode=twoCamMode, deRotate=doDerotate, Ibias=0)
#p.makeDiffIms(imRange=imRange, showPlots=False, maskfile='manualMasksNoc2016.npz')
#p.makeDiffIms(showPlots=False)
# mat1 = np.identity(4)
# mat2 = np.identity(4)
# p.makeDiffIms(imRange=imRange, showPlots=False, calMat_PreRot=mat1, calMat_PostRot=mat2)

# p.makeDiffIms(imRange=imRange, showPlots=False, twoCamMode=twoCamMode, deRotate=doDerotate, doFlats=True)

I = p.StokesImsSummed[:,:, 0]
Q = p.StokesImsSummed[:,:, 1]
U = p.StokesImsSummed[:,:, 2]
pQ = Q / I
pU = U / I
p2_frac = np.sqrt(pQ**2 + pU**2)
p_frac = np.sqrt(Q**2 + U**2) / I
P = np.sqrt(Q**2 + U**2)

#p.plotSingleIm(imType='p', clim=[-0.05, 0.4])
#p.plotSingleIm(imType='p', cmap='gray', crop=0.5)

im = p2_frac
# im = np.clip(im,-0.1,0.4)
# # im = U
# #im = p.polzRatioImsSummed
curHDU = fits.PrimaryHDU()
curHDU.data = im
curHDU.writeto('out.fits', overwrite=True)

curHDU = fits.PrimaryHDU()
curHDU.data = np.transpose(p.StokesImsSummed)
curHDU.writeto('outStokes.fits', overwrite=True)


# plt.imshow(im)
# plt.imshow(im, cmap='gray')
p.plotSingleIm(imType='p',crop=0.5, cmap='viridis')
plt.pause(0.001)


# # Show polarisation ratio images
# im_f = p.polzRatioImsSummed
# plt.imshow(im_f)
# plt.imshow(im_f[48:208, 48:208])
# #plt.pause(0.001)

# fits.writeto('out.fits', data, header)

"""
dataPath = '/import/silo4/snert/VAMPIRESData_201609/Analysis/MWC758/'
filePref = 'cube_MWC758_01-02Combined_20160919_750-50_EmptySlot_'
nSubFiles = 172 *4 # NB actual files, not filenums. So 4x num raw files

import simplePDIFull_01 as s
cubeInfoFile = 'cubeinfoOct2016.idlvar'
p = s.pdiImages()
p.processRawFiles(dataPath, filePref, cubeInfoFile, nSubFiles,saveFilePref =
    'allSummedImsCube_maxMeth_', method='max')
"""

print('Done.')



"""
from rebin import rebin
rebinSize = 4

pQ=p.plotSingleIm(imType='pQ',crop=0.5, cmap='viridis')
pU=p.plotSingleIm(imType='pU',crop=0.5, cmap='viridis')

pQ_r = rebin(pQ, rebinSize)
pU_r = rebin(pU, rebinSize)
p_r = np.sqrt(pQ_r**2 + pU_r**2)
psi_r = 0.5*np.arctan(pU_r/pQ_r)
psi = 0.5*np.arctan(pU/pQ)

psi_r = -psi_r + np.pi/2 #0 should be up, not right...?

plt.clf()
plt.imshow(p_r)
# Use quiver but give the angles explicitly, and magnitudes by setting one input to fractional P
# plt.quiver(p_r, 0, angles=psi_r/np.pi*180, pivot='mid', headwidth=0, scale=0.4)
plt.quiver(p_r, 0, angles=psi_r/np.pi*180, pivot='mid', headwidth=0)
"""



