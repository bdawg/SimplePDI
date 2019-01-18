"""
Script to run simplePDI with various tasks
"""
from simplePDI import *



######################################################################################################
# Usage for reading *raw* files (no IDL pre-cubing)
IDLPreCubed = False

dataPath = "/Users/bnorris/DontBackup/vampdata_Oct2018/"
filePref = 'ABAur_01_20181018_750-50_EmptySlot_0'
nSubFiles = 40 *2*2
saveFilePref = 'allSummedImsCube_MAXmeth_darksubt_DONEWITHSIMPLEPDI_'
nSubFiles = 4 *2*2

dataPath = "/Users/bnorris/DontBackup/pcal_20180712/"
filePref = 'pcal20180712_set2_nopol_750-50_EmptySlot_'
nSubFiles = 1 *2*2
saveFilePref = 'tempsave'
darkFilename = '../SimplePDI_DATA/summedDarks_pcal20180712_darks_775-50_Mirror_.npz'

# dataPath = "/Users/bnorris/DontBackup/vampdata_polzEngStars/"
# filePref = 'RhoBoo_02_20180625_750-50_EmptySlot_0'
# nSubFiles = 4 *2*2
# saveFilePref = 'tempsave'
# darkFilename = '../SimplePDI_DATA/summedDarks_dark_EM10_20180622_625-50_Mirror_0.npz'
# # From /Users/bnorris/DontBackup/vampdata_polzEngStars/dark_EM10_20180622_625-50_Mirror_...


startFileNum = 0
centroidThresh = 5000 #Values below this clipped to zero for centroiding

# What to do about bg subtraction?
# 0 - don't subtract anything
# 1 - subtract specified bias values, as per bgVals tuple
# 2 - subtract actual darkframes from file specified in darkFilename
# 3 - use a sky annulus, as per skyAnnInner and Outer values. (Careful - halo is big!)
bgMode = 2
bgVals = (176.5, 178.9) # (chan1, chan2)

# Be careful with outer radius to avoid satellite spots!
skyAnnInner = 90
skyAnnOuter = 120


# darkFilename = None
# bgMode = 0




######################################################################################################
# Settings for making the differential images

twoCamMode = True # If true, just reflect chan2 horizontally, not rotate.
doDerotate = True
# bsRot = 89.93

# dataPath = '../SimplePDI_DATA/'
# dataPath = '/Users/bnorris/DontBackup/vampdata_Oct2018/'
#
#
# loadFilename = 'allSummedImsCube_MAXmeth_darksubt_ABAur_01_20181017_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_ABAur201612_52acqs_COMmethABAur_02_20161215_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_MAXmeth_darksubt_ABAur_01_20181018_750-50_EmptySlot_0.npz'




######################################################################################################
# Do it



p = pdiImages(outDir='../SimplePDI_DATA/')


p.processRawFiles(dataPath, filePref, nSubFiles, startFileNum=startFileNum,
                  saveFilePref=saveFilePref, centroidThresh=centroidThresh,
                  method='com', twoCamMode = True, bgMode=bgMode, darkFilename=darkFilename,
                  showAllIms=False, previewCropRad=25, comRad=4,
                  luckyCriteria=None, luckyPercent=99, singleFileOnly=True)

# p.loadCube(dataPath, loadFilename, twoCamMode = twoCamMode)

p.makeDiffIms(showPlots=False, twoCamMode=twoCamMode, deRotate=doDerotate)

# p.plotSingleIm(imType='p',crop=0.6, cmap='viridis')
#
# p.plotStokesIms(crop=0.6, saveFITS=True) #, savePath=dataPath)

p.partialDiffcal('1a', showImageCrop=0.5, showImageType='pQ', clim=[0,1])

im=p.plotSingleIm('partial', imType='pQ', clim=[0,1])
p.measureAperture(im, 20, showIm=True, mean=True, showImSize=48)

################################## TESTING ###################################

# # Get the first HWP set images:
# cur=p.allPolzstateIms[0]
# im_c0 = cur.h0c0l0
# im_c1 = cur.h0c1l0
# im_q = im_c0 - im_c1
#
# plt.figure(3)
# plt.imshow(im_q)
# plt.colorbar()