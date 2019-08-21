from preprocVampiresFunc import *

dataPath = '/import/silo4/snert/VAMPIRESData/201906/20190624/'
outputPath = '/import/silo4/snert/VAMPIRESData/201906/201906_Yang_preproc/'

ALL_filePrefs = []
ALL_nSubFiles = []
ALL_outFilePrefs = []


# Each list entry is a 3-elementt list of [filePrefs, nSubFilesAll, outFilePref]
allSettings = [
    # [ [''], [], [''] ],
    [ ['HD141569_01_20190625_750-50_EmptySlot_0', 'HD141569_01__RESTARTED___20190625_750-50_EmptySlot_0', 
    	'HD141569_02_20190625_750-50_EmptySlot_0', 'HD141569_04__RESTARTED___20190625_750-50_EmptySlot_0'], 
    	[52, 28, 72, 76], ['HD141569_20190625_750-50_Open_'] ],
    [ ['HIP77410_01_20190625_750-50_EmptySlot_0'], [24], ['HIP77410_20190625_750-50_Open_'] ],
    [ ['HD149914_01_20190625_750-50_EmptySlot_0', 'HD149914_01__RESTARTED___20190625_750-50_EmptySlot_0'], [68, 24], ['HD149914_20190625_750-50_Open_'] ],
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


















