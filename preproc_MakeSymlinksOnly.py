"""
Pre-process VAMPIRES data:
Extracts interleaved polarisation states, discards leading frames, fixes filenames and
saves to new FITS files.

The output files will be of format outFilePref_[serial number]_[camera]_[FLC].fits
[camera] is labelled as '1' or '2' and
[FLC] is labelled as 'A' or 'B'.
"""

# Lists of multiple filePrefs and nSubFiles may be used, and sets will be concatenated.
# NB separate subfiles for each state, not original filenums.

# dataPath = '/import/silo4/snert/VAMPIRESData/201902/20190226/'
# filePrefs = ['HD34700_06_20190226_750-50_EmptySlot_0',
#             'HD34700_08_20190226_750-50_EmptySlot_0',
# 			'HD34700_08__RESTARTED___20190226_750-50_EmptySlot_0',
# 			'HD34700_08__RESTARTED____RESTARTED___20190226_750-50_EmptySlot_0']
# nSubFilesAll = [8,4,4,4]
# outputPath = '/import/silo4/snert/VAMPIRESData/201902/20190226/linked/'
# outFilePref = 'HD34700_06-08_20190226_750-50_EmptySlot_0'

# dataPath = '/Volumes/BigVampData/201902/20190225//'
# filePrefs = ['',
#             '',
# 			'',
# 			'']
# nSubFilesAll = [,,,]
# outputPath = '/Volumes/BigVampData/201902/20190225/linked/'
# outFilePref = 'HD34700_06-08_20190226_750-50_EmptySlot_0'

dataPath = '/Volumes/BigVampData/SELECTED_COPIES/from20170116/'
filePrefs = ['ABAur_01_20170116_750-50_EmptySlot_0',
            'ABAur_02_20170116_750-50_EmptySlot_0',
			'ABAur_03_20170116_750-50_EmptySlot_0']
nSubFilesAll = [68,176,88]
outputPath = '/Volumes/BigVampData/SELECTED_COPIES/from20170116/'
outFilePref = 'ABAur_01-03_20170116_750-50_EmptySlot_0'

dataPath = '/Volumes/BigVampData/201801/20180108/'
filePrefs = ['ABAur_01_20180108_750-50_EmptySlot_0',
            'ABAur_02_20180108_750-50_EmptySlot_0']
nSubFilesAll = [108,84]
outputPath = '/Volumes/BigVampData/201801/20180108/linked/'
outFilePref = 'ABAur_01-02_20180108_750-50_EmptySlot_0'

dataPath = '/Volumes/BigVampData/201810/20181022/'
filePrefs = ['ABAur_02_20181022_750-50_EmptySlot_0',
            'ABAur_03_20181022_750-50_EmptySlot_0',
            'ABAur_03__RESTARTED___20181022_750-50_EmptySlot_0' ]
nSubFilesAll = [24,24,20]
outputPath = '/Volumes/BigVampData/201810/20181022/linked/'
outFilePref = 'ABAur_02-03_combined_20181022_750-50_EmptySlot_0'

dataPath = '/Volumes/BigVampData/201810/20181018/'
filePrefs = ['ABAur_01_20181018_750-50_EmptySlot_0',
            'ABAur_03_20181018_750-50_EmptySlot_0',
            ]
nSubFilesAll = [40,4]
outputPath = '/Volumes/BigVampData/201810/20181018/linked/'
outFilePref = 'ABAur_01-03_goodfiles_20181018_750-50_EmptySlot_0'

dataPath = '/Volumes/BigVampData/201812/20181214/'
filePrefs = ['ABAur_05_20181215_750-50_EmptySlot_0',
            'ABAur_05__RESTARTED___20181215_750-50_EmptySlot_0',
            'ABAur_06_20181215_750-50_EmptySlot_0',
            'ABAur_07_20181215_750-50_EmptySlot_0',
            'ABAur_07__RESTARTED___20181215_750-50_EmptySlot_0']
nSubFilesAll = [8,76,160,72,120]
outputPath = '/Volumes/BigVampData/201812/20181214/linked/'
outFilePref = 'ABAur_05-07_20181214_750-50_EmptySlot_0'

dataPath = '/Volumes/BigVampData/201810/20181023/'
filePrefs = ['mwc758_01_20181023_Open_EmptySlot_0',
            'mwc758_02_20181023_Open_EmptySlot_0',
            'mwc758_02__RESTARTED___20181023_Open_EmptySlot_0',
            'mwc758_02__RESTARTED____RESTARTED___20181023_Open_EmptySlot_0']
nSubFilesAll = [92,64,96,60]
outputPath = '/Volumes/BigVampData/201810/20181023/linked/'
outFilePref = 'mwc758_01-02_20181023_Open_EmptySlot_0'





dryRun = False
fileExtn = '.fits'
startFileNum = 0

##############################################################################################
import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
from copy import copy
import time
import os


# Indexes are [:, :, hwpSet, HWP, Channel (camera), FLCstate]
# So each set of 4 VAMPIRES on-sky files corresponds to 1 hwpSet.
curOutFnum = 0
for fSet in range(len(filePrefs)):
    filePref = filePrefs[fSet]
    nSubFiles = nSubFilesAll[fSet]

    for f in range(0, nSubFiles):
        curFileNum = f

        for c in range(0, 2):
            curChan = c

            # Generate current filename
            if curChan == 0:
                curCamStr = '_cam1'
                #curCamStrOut = '_1'
            else:
                curCamStr = '_cam2'
                #curCamStrOut = '_2'

            curFilenumStr = '%d' % curFileNum
            curFilename = dataPath + filePref + curFilenumStr + curCamStr + fileExtn
            print('From file %s' % curFilename)
            outFilename = outputPath + outFilePref + str(curOutFnum) + curCamStr + '.fits'
            print('To file %s' % outFilename)

            cmdStr = 'ln -s ' + curFilename + ' ' + outFilename
            print(cmdStr)
            print(' ')
            if not dryRun:
                os.system(cmdStr)

        curOutFnum = curOutFnum + 1

























