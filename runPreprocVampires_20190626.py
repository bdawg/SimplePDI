from preprocVampiresFunc import *

dataPath = '/import/silo4/snert/VAMPIRESData/201906/20190626/'
outputPath = '/import/silo4/snert/VAMPIRESData/201906/201906_Yang_preproc/'

ALL_filePrefs = []
ALL_nSubFiles = []
ALL_outFilePrefs = []


# Each list entry is a 3-elementt list of [filePrefs, nSubFilesAll, outFilePref]
allSettings = [
    # [ [''], [], [''] ],
    [ ['HD145718_01_20190626_750-50_EmptySlot_0', 'HD145718_02_20190626_750-50_EmptySlot_0'], [4, 52], ['HD145718_20190626_750-50_Open_'] ],
    [ ['HIP79366_01_20190626_750-50_EmptySlot_0'], [8], ['HIP79366_20190626_750-50_Open_'] ],
    [ ['HD139614_01_20190626_750-50_EmptySlot_0', 'HD139614_01__RESTARTED___20190626_750-50_EmptySlot_0'], [4, 80], ['HD139614_20190626_750-50_Open_'] ],
    [ ['HD140243_01_20190626_750-50_EmptySlot_0'], [12], ['HD140243_01_20190626_750-50_Open_'] ],
    [ ['V4046Sgr_01_20190626_750-50_EmptySlot_0', 'V4046Sgr_01__RESTARTED___20190626_750-50_EmptySlot_0'], [52, 20], ['V4046Sgr_01_20190626_750-50_Open_'] ],
    [ ['HIP89105_01_20190626_750-50_EmptySlot_0'], [12], ['HIP89105_01_20190626_750-50_Open_'] ],
]



# allSettings = [
#     [['V1247Ori_01_20190320_Open_EmptySlot_0'], [2], ['V1247Ori_01_20190320_Open_']],
#     [['HD37845_01_20190320_Open_EmptySlot_0'], [4], ['HD37845_01_20190320_Open_']],
#     [['GWOri_01_20190320_Open_EmptySlot_0', 'GWOri_03_20190320_Open_EmptySlot_0'], [2, 4], ['GWOri_20190320_Open_']],
#     ]


from preprocVampiresFunc import *

for curSettings in allSettings:
    filePrefs = curSettings[0]
    nSubFilesAll = curSettings[1]
    outFilePref = curSettings[2][0]

    preprocVampiresFunc(dataPath, filePrefs, nSubFilesAll, outputPath, outFilePref)


print('Finished!')


















