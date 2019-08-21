from preprocVampiresFunc import *

dataPath = '/import/silo4/snert/VAMPIRESData/201903/20190320/'
outputPath = '/import/silo4/snert/VAMPIRESData/201903/20190320_preproc/'

ALL_filePrefs = []
ALL_nSubFiles = []
ALL_outFilePrefs = []


# Each list entry is a 3-elementt list of [filePrefs, nSubFilesAll, outFilePref]
allSettings = [
    # [ [''], [], [''] ],
    [ ['domeflat_750_em300_1s_20190320_750-50_EmptySlot_0'], [4], ['domeflat_750_em300_1s_20190320_750-50_'] ],
    [ ['domeflat_750_em1000_20ms_20190320_750-50_EmptySlot_0'], [16], ['domeflat_750_em1000_20ms_20190320_750-50_']],
    [ ['ABAur_01_20190320_750-50_EmptySlot_0'], [104], ['ABAur_01_20190320_750-50_'] ],
    [ ['HD31706_01_20190320_750-50_EmptySlot_0'], [8], ['HD31706_01_20190320_750-50_'] ],
    [ ['HD31706_02_20190320_750-50_EmptySlot_0'], [8], ['HD31706_02_20190320_750-50_'] ],
    [['V1247Ori_01_20190320_Open_EmptySlot_0'], [20], ['V1247Ori_01_20190320_Open_']],
    [['HD37845_01_20190320_Open_EmptySlot_0'], [4], ['HD37845_01_20190320_Open_']],
    [['GWOri_01_20190320_Open_EmptySlot_0', 'GWOri_03_20190320_Open_EmptySlot_0'], [4, 20], ['GWOri_20190320_Open_']],
    [['HD244992_01_20190320_Open_EmptySlot_0'], [8], ['HD244992_01_20190320_Open_']],
    [['TWHya_01_20190320_Open_EmptySlot_0', 'TWHya_01__RESTARTED___20190320_Open_EmptySlot_0'], [14, 20], ['TWHya_20190320_Open_']],
    [['CD337042_01_20190320_Open_EmptySlot_0'], [8], ['CD337042_01_20190320_Open_']],
    [['TWA7_02_20190320_Open_EmptySlot_0', 'TWA7_02__RESTARTED___20190320_Open_EmptySlot_0'], [20, 12], ['TWA7_20190320_Open_']],
    [['HD98800_01_20190320_Open_EmptySlot_0'], [36], ['HD98800_01_20190320_Open_']],
    [['CD258554_01__RESTARTED___20190320_Open_EmptySlot_0'], [4], ['CD258554_20190320_Open_EmptySlot_0']],
    [['SAO206462_01_20190320_Open_EmptySlot_0', 'SAO206462_01__RESTARTED___20190320_Open_EmptySlot_0'], [20, 12], ['SAO206462_20190320_Open_']],
    [['HD135562_01_20190320_Open_EmptySlot_0'], [8], ['HD135562_01_20190320_Open_']],
    [['M5_6_20190320_Open_EmptySlot_0'], [4], ['M5_6_20190320_Open_']],
    [['M5_7_20190320_Open_EmptySlot_0'], [4], ['M5_7_20190320_Open_']],
    [['HD154445_01_20190320_Open_EmptySlot_0'], [8], ['HD154445_01_20190320_Open_']],
    [['HD154445_02_20190320_750-50_EmptySlot_0'], [8], ['HD154445_02_20190320_750-50_']],
    [['HD154445_03_20190320_750-50_EmptySlot_0'], [8], ['HD154445_03_20190320_750-50_']],
    [['HD143006_01_20190320_Open_EmptySlot_0'], [36], ['HD143006_01_20190320_Open_']],
    [['BD224047_01_20190320_Open_EmptySlot_0'], [8], ['BD224047_01_20190320_Open_']],
    [['twilightFlat_02_20190320_Open_EmptySlot_0'], [4], ['twilightFlat_02_20190320_Open_']],
    [['dark_1000ms_em300_20190320_Open_Mirror_0'], [4], ['dark_1000ms_em300_20190320_']],
    [['dark_100ms_em300_20190320_Open_Mirror_0'], [20], ['dark_100ms_em300_20190320_']],
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


















