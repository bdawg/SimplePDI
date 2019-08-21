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
centroidThresh = 5000 #Values below this clipped to zero for centroiding

# dataPath = "/Users/bnorris/DontBackup/vampdata_polzEngStars/"
# filePref = 'RhoBoo_02_20180625_750-50_EmptySlot_0'
# nSubFiles = 4 *2*2
# saveFilePref = 'tempsave'
# darkFilename = '../SimplePDI_DATA/summedDarks_dark_EM10_20180624_Open_Mirror_0.npz'
# # From /Users/bnorris/DontBackup/vampdata_polzEngStars/dark_EM10_20180622_625-50_Mirror_...
# centroidThresh = 5000 #Values below this clipped to zero for centroiding
#
# dataPath = "/Users/bnorris/DontBackup/vampdata_polzEngStars/"
# filePref = 'HR8799_05_nicePSF_20170912_750-50_EmptySlot_0'
# nSubFiles = 4 *2*2
# saveFilePref = 'tempsave'
# darkFilename = '../SimplePDI_DATA/summedDarks_dark_EM10_20180624_Open_Mirror_0.npz'
# # darkFilename = '../SimplePDI_DATA/summedDarks_darks_10ms_em20_01_20170903_750-50_Mirror_0.npz'
# centroidThresh = 500 #Values below this clipped to zero for centroiding
#
# dataPath = "/Users/bnorris/DontBackup/vampdata_polzEngStars/"
# filePref = '63Cet_01_20170912_750-50_EmptySlot_0'
# nSubFiles = 16 *2*2
# saveFilePref = 'tempsave'
# darkFilename = '../SimplePDI_DATA/summedDarks_dark_EM10_20180624_Open_Mirror_0.npz'
# # From /Users/bnorris/DontBackup/vampdata_polzEngStars/dark_EM10_20180622_625-50_Mirror_...
# centroidThresh = 500 #Values below this clipped to zero for centroiding

dataPath = "/Users/bnorris/DontBackup/vampdata_polzEngStars/"
filePref = 'altair_02_20180624_750-50_EmptySlot_0'
nSubFiles = 4 *2*2
saveFilePref = 'tempsave_'
darkFilename = '../SimplePDI_DATA/summedDarks_dark_EM00_20180624_Open_Mirror_0.npz'
centroidThresh = 500 #Values below this clipped to zero for centroiding

dataPath = "/Users/bnorris/DontBackup/vampdata_polzEngStars/"
filePref = 'GJ504_01_20170619_750-50_EmptySlot_0'
nSubFiles = 8 *2*2
saveFilePref = 'tempsave_'
darkFilename = '../SimplePDI_DATA/summedDarks_darks_10ms_em20_01_20170903_750-50_Mirror_0.npz'
centroidThresh = 500 #Values below this clipped to zero for centroiding

dataPath = "/Users/bnorris/DontBackup/vampdata_polzEngStars/"
filePref = 'phiAqr_750_20180622_750-50_EmptySlot_0'
nSubFiles = 4 *2*2
saveFilePref = 'tempsave_'
darkFilename = '../SimplePDI_DATA/summedDarks_dark_EM10_20180622_625-50_Mirror_0.npz'
centroidThresh = 500 #Values below this clipped to zero for centroiding

dataPath = "/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201704/20170418/"
filePref = 'HD183143_03_20170418_750-50_EmptySlot_0'
nSubFiles = 16 *2*2
saveFilePref = 'tempsave_'
darkFilename = '../SimplePDI_DATA/summedDarks_darks_10ms_em20_01_20170903_750-50_Mirror_0.npz'
centroidThresh = 500 #Values below this clipped to zero for centroiding

# dataPath = "/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201709/20170903/"
# filePref = 'HD183143_01_20170903_750-50_EmptySlot_0'
# nSubFiles = 92 *2*2
# saveFilePref = 'tempsave_'
# darkFilename = '../SimplePDI_DATA/summedDarks_darks_10ms_em300_01_20170903_750-50_Mirror_0_16Files.npz'
# centroidThresh = 1000 #Values below this clipped to zero for centroiding
#
# dataPath = "/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201709/20170904/"
# filePref = 'HD183143_01_20170904_750-50_EmptySlot_0'
# nSubFiles = 72 *2*2
# saveFilePref = 'tempsave_'
# darkFilename = '../SimplePDI_DATA/summedDarks_darks_10ms_em300_01_20170903_750-50_Mirror_0_16Files.npz'
# centroidThresh = 1000 #Values below this clipped to zero for centroiding
#
# dataPath = "/Volumes/silo4/snert/VAMPIRESData/201806/20180622/"
# filePref = 'phiAqr_750_20180622_750-50_EmptySlot_0'
# nSubFiles = 4 *2*2
# saveFilePref = 'allImsCube_50pc_'
# darkFilename = '../SimplePDI_DATA/summedDarks_dark_EM10_20180622_625-50_Mirror_0.npz'
# centroidThresh = 1000 #Values below this clipped to zero for centroiding
#
# dataPath = "/Volumes/silo4/snert/VAMPIRESData/201812/20181214/"
# filePref = 'HD282411_03_20181215_750-50_EmptySlot_0'
# nSubFiles = 8 *2*2
# saveFilePref = 'allImsCube_50pc_'
# darkFilename = '../SimplePDI_DATA/summedDarks_dark_em300_20180628_625-50_Mirror_0.npz'
# centroidThresh = 1000 #Values below this clipped to zero for centroiding

dataPath = "/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201705/20170531/"
filePref = 'vega_02_20170531_750-50_18holeNudged_0'
nSubFiles = 24 *2*2
saveFilePref = 'tempsave_'
darkFilename = '../SimplePDI_DATA/summedDarks_darks_10ms_em20_01_20170903_750-50_Mirror_0.npz'
centroidThresh = 500 #Values below this clipped to zero for centroiding

dataPath = "/Volumes/RedDisk/VAMPIRES_AllData_2018/201812/20181215/"
filePref = 'HD44612_01_20181215_750-50_EmptySlot_0'
nSubFiles = 16 *2*2
saveFilePref = 'tempsave_'
darkFilename = '../SimplePDI_DATA/summedDarks_SKIES_02_20181215_750-50_EmptySlot_0.npz'
centroidThresh = 100 #Values below this clipped to zero for centroiding

dataPath = "/Volumes/RedDisk/VAMPIRES_AllData_2018/201810/20181022/"
filePref = 'HD30454_03_20181022_750-50_EmptySlot_0'
nSubFiles = 40 *2*2
saveFilePref = 'allImsCube_50pc_'
darkFilename = '../SimplePDI_DATA/summedDarks_darks_5ms_128px_em25_20181022_750-50_Mirror_0.npz'
centroidThresh = 500 #Values below this clipped to zero for centroiding





dataPath = "/Volumes/BigVampData/201902/20190225/linked/"
filePref = 'MWC480_02-04_20190225_broadband_EmptySlot_0'
nSubFiles = 268 *2*2
saveFilePref = 'allImsCube_'
darkFilename = '../SimplePDI_DATA/summedDarks_dark_256_10_em300_02_20190225_Open_Mirror_0.npz'
centroidThresh = 500 #Values below this clipped to zero for centroiding

dataPath = "/Volumes/BigVampData/201903/20190320/"
filePref = 'ABAur_01_20190320_750-50_EmptySlot_0'
nSubFiles = 104 *2*2
saveFilePref = 'allImsCube_'
darkFilename = '../SimplePDI_DATA/summedDarks_dark_20ms_em300_20190320_Open_Mirror_0.npz'
centroidThresh = 1000 #Values below this clipped to zero for centroiding

dataPath = "/Volumes/BigVampData/201903/20190320/"
filePref = 'domeflat_750_em1000_20ms_20190320_750-50_EmptySlot_0'
nSubFiles = 16 *2*2
saveFilePref = 'allImsCube_'
darkFilename = '../SimplePDI_DATA/summedDarks_dark_20ms_em300_512px_20190320_Open_Mirror_0.npz'
centroidThresh = 1000 #Values below this clipped to zero for centroiding

dataPath = "/Volumes/BigVampData/201903/20190320/"
filePref = 'domeflat_750_em300_1s_20190320_750-50_EmptySlot_0'
nSubFiles = 4 *2*2
saveFilePref = 'allImsCube_'
darkFilename = '../SimplePDI_DATA/summedDarks_dark_1000ms_em300_20190320_Open_Mirror_0.npz'
centroidThresh = 1000 #Values below this clipped to zero for centroiding

# start image _089_cam2 at 5:45pm
# 97_cam1 att 5:52

startFileNum = 0

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

luckyCriteria = None
# luckyCriteria = 'l2flux'
luckyPercent = 95


# darkFilename = None
# bgMode = 0


showAllIms = False
previewCropRad = 250 #25
comRad=4
comRad = 8 # Usually 4 is good, but more if saturated...

######################################################################################################
# Settings for making the differential images

twoCamMode = True # If true, just reflect chan2 horizontally, not rotate.
doDerotate = True
# bsRot = 89.93

# dataPah = '../SimplePDI_DATA/'
# dataPath = '/Users/bnorris/DontBackup/vampdata_Oct2018/'
#
#
# loadFilename = 'allSummedImsCube_MAXmeth_darksubt_ABAur_01_20181017_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_ABAur201612_52acqs_COMmethABAur_02_20161215_750-50_EmptySlot_0.npz'
# loadFilename = 'allSummedImsCube_MAXmeth_darksubt_ABAur_01_20181018_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_pcal20180712_set2_nopol_750-50_EmptySlot_.npz'
# loadFilename = 'allImsCube_50pc_altair_02_20180624_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_63Cet_01_20170912_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_GJ504_01_20170619_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_HR8799_05_nicePSF_20170912_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_lucky50pc_RhoBoo_02_20180625_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_HD155528_02_20180626_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_HD154445_02_20180626_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_phiAqr_750_20180622_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_HD206893_02_20180623_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_HD282411_03_20181215_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_HD43384_02_20181215_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_HD43384_02_20181023_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_HD183143_01_20170903_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_HD183143_01_20170904_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_HD183143_02_20180622_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_HD183143_03_20170418_750-50_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_HD198478_04_20180623_750-50_EmptySlot_0.npz'

loadFilename = 'allImsCube_MWC480_02-04_20190225_broadband_EmptySlot_0.npz'
# loadFilename = 'allImsCube_50pc_HD44612_01_20181215_750-50_EmptySlot_0.npz'
loadFilename = 'allImsCube_ABAur_01-02_20180108_750-50_EmptySlot_0.npz' # The weird circle-smeared one

######################################################################################################
# Do it



p = pdiImages(outDir='../SimplePDI_DATA/')


p.processRawFiles(dataPath, filePref, nSubFiles, startFileNum=startFileNum,
                  saveFilePref=saveFilePref, centroidThresh=centroidThresh,
                  method='none', twoCamMode = True, bgMode=bgMode, darkFilename=darkFilename,
                  showAllIms=showAllIms, previewCropRad=previewCropRad, comRad=comRad,
                  luckyCriteria=luckyCriteria, luckyPercent=luckyPercent, singleFileOnly=False)

# p.loadCube(dataPath, loadFilename, twoCamMode = twoCamMode)

p.makeDiffIms(showPlots=False, twoCamMode=twoCamMode, deRotate=doDerotate, darkOffset=None)

p.plotStokesIms(crop=0.25, apertureRad=20)

# p.doMultiCals(apertureRad=20)
# p.doMultiCals(apertureRad=8)



# p.makeDiffIms(showPlots=False, twoCamMode=twoCamMode, deRotate=doDerotate, darkOffset=None, paOffset = 54.)
# p.plotStokesIms(crop=0.25, apertureRad=20)


# For lab data 20180712:
# p.processRawFiles(dataPath, filePref, nSubFiles, startFileNum=startFileNum,
#                   saveFilePref=saveFilePref, centroidThresh=centroidThresh,
#                   method='com', twoCamMode = True, bgMode=bgMode, darkFilename=darkFilename,
#                   showAllIms=False, previewCropRad=25, comRad=4,
#                   luckyCriteria=None, luckyPercent=99, singleFileOnly=True)

# p.loadCube(dataPath, loadFilename, twoCamMode = twoCamMode)

# p.makeDiffIms(showPlots=False, twoCamMode=twoCamMode, deRotate=doDerotate)
#
# # p.plotSingleIm(imType='p',crop=0.6, cmap='viridis')
# #
# # p.plotStokesIms(crop=0.6, saveFITS=True) #, savePath=dataPath)
#
# p.partialDiffcal('1a', showImageCrop=0.5, showImageType='pQ', clim=[0,1])
#
# im=p.plotSingleIm('partial', imType='pQ', clim=[0,1])
# p.measureAperture(im, 16, showIm=True, mean=True, showImSize=48)

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