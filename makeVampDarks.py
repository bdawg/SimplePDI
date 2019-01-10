import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt
plt.ion()

fileExtn = '.fits'
dataPath = '/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201709/20170903/'
filePref = 'darks_10ms_em300_01_20170903_750-50_Mirror_0'
nSubFiles = 32 # e.g. 8 acquisitions with 2 cameras = 16 subfiles
startFileNum = 0
HAFileformat = False

dataPath = '/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201701/20170119/'
filePref = 'darks_20170118_512_em1000_18000us_20170119_750-50_Mirror_0'
nSubFiles = 16 *2 # e.g. 8 acquisitions with 2 cameras = 16 subfiles


dataPath = '/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201801/20180108/'
dataPath = '/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201801/20180109/'
filePref = 'dark_PupilWheelBlock_05s300_01_HaDifferential-'
filePref = 'dark_PupilWheelBlock_05s300_512_02_HaDifferential-'
filePref = 'dark_CameraShutter_05s300_256_HaDifferential-'
filePref = 'dark_CameraShutter_05s300_512_HaDifferential-'
fileSuf = '_Open_AnnulusNudged_'
fileSuf = '_Open_EmptySlot_'
startFileNum = 0
nSets = 4 # 1 Set is two states x two cams
HAFileformat = True

dataPath = '/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201801/tempcopy/'
filePref = 'dark_PupilWheelBlock_05s300_256_longset_COMBINED__20180109_Open_EmptySlot_0'
nSubFiles = 100 *2 # e.g. 8 acquisitions with 2 cameras = 16 subfiles
HAFileformat = False

dataPath = '/Users/bnorris/DontBackup/vampdata_Oct2018/'
filePref = 'darks_40ms_em300_20181017_750-50_Mirror_0'
nSubFiles = 8 *2 # e.g. 8 acquisitions with 2 cameras = 16 subfiles
HAFileformat = False


useMedian = True # Take the median of all files, not the mean
showSEMMap = True
saveData = True
saveFilePref = 'summedDarks_'

plt.figure()
if HAFileformat:
    # Filter arrangement, listed as [cam1, cam2]
    state1 = ('cont', 'Ha')
    state2 = ('Ha', 'cont')

    # Read a FITS file to get sizes
    curFilenumStr = '%d' % startFileNum
    curCamStr = '_cam1'
    curStateStr = 'state1'
    curFilename = dataPath + filePref + curStateStr + fileSuf + curFilenumStr + curCamStr + fileExtn
    hdulist = fits.open(curFilename)
    curHDU = hdulist[0]
    curCube = np.transpose(curHDU.data)
    nFrms = curCube.shape[2]
    dim = curCube.shape[0]

    # Indexes are [:, :, Set+State, Channel (camera)]
    allSummedIms = np.zeros([dim, dim, nSets*2, 2])

    curSetState = 0
    for f in range(0, nSets):
        curFileNum = f
        curFilenumStr = '%d' % curFileNum

        for s in range(0, 2):
            curState = s
            if curState == 0:
                curStateStr = 'state1'
            else:
                curStateStr = 'state2'

            for c in range(0, 2):
                curChan = c
                if curChan == 0:
                    curCamStr = '_cam1'
                else:
                    curCamStr = '_cam2'

                curFilename = dataPath + filePref + curStateStr + fileSuf + curFilenumStr \
                              + curCamStr + fileExtn
                print('Reading file %s' % curFilename)
                hdulist = fits.open(curFilename)
                curHDU = hdulist[0]
                curCube = np.transpose(curHDU.data)
                goodframes = curCube[:, :, 2:nFrms]  # Discard 1st 2 frames
                curDark = np.mean(goodframes, axis=2)
                allSummedIms[:, :, curSetState, curChan] = curDark

                plt.clf()
                plt.imshow(curDark)
                plt.colorbar()
                plt.pause(0.001)

            curSetState = curSetState + 1


else:
    # Read a FITS file to get sizes
    curFilenumStr = '%d' % startFileNum
    curCamStr = '_cam1'
    curFilename = dataPath + filePref + curFilenumStr + curCamStr + fileExtn
    hdulist = fits.open(curFilename)
    curHDU = hdulist[0]
    curCube = np.transpose(curHDU.data)
    nFrms = curCube.shape[2]
    dim = curCube.shape[0]

    allSummedIms = np.zeros([dim, dim, nSubFiles // 2, 2])
    curSet = 0
    for f in range(0, nSubFiles // 2):
        curFileNum = f

        for c in range(0, 2):
            curChan = c

            # Generate current filename
            if curChan == 0:
                curCamStr = '_cam1'
            else:
                curCamStr = '_cam2'

            curFilenumStr = '%d' % curFileNum
            curFilename = dataPath + filePref + curFilenumStr + curCamStr + fileExtn

            print('Reading file %s' % curFilename)
            hdulist = fits.open(curFilename)
            curHDU = hdulist[0]
            curSuperCube = np.transpose(curHDU.data)
            goodframes = curSuperCube[:, :, 2:nFrms] # Discard 1st 2 frames
            curDark = np.mean(goodframes, axis=2)
            allSummedIms[:, :, curSet, curChan] = curDark

            if curChan == 1:
                curSet = curSet + 1


if useMedian:
    finalDarks = np.median(allSummedIms, axis=2)
else:
    finalDarks = np.mean(allSummedIms, axis=2)
sigmaDarks = np.std(allSummedIms, axis=2) / np.sqrt(nFrms-1)

if showSEMMap:
    plt.figure()
    plt.subplot(2, 2, 1)
    plt.imshow(finalDarks[:, :, 0])
    plt.colorbar()
    plt.title('Camera 1 Dark')
    plt.subplot(2, 2, 2)
    plt.imshow(finalDarks[:, :, 1])
    plt.colorbar()
    plt.title('Camera 2 Dark')

    plt.subplot(2, 2, 3)
    plt.imshow(sigmaDarks[:, :, 0])
    plt.colorbar()
    plt.title('Camera 1 SEM map')
    plt.subplot(2, 2, 4)
    plt.imshow(sigmaDarks[:, :, 1])
    plt.colorbar()
    plt.title('Camera 2 SEM map')

else:
    plt.figure(figsize=(18,6))
    plt.subplot(1, 2, 1)
    plt.imshow(finalDarks[:, :, 0])
    plt.colorbar()
    plt.title('Camera 1 Dark')
    plt.subplot(1, 2, 2)
    plt.imshow(finalDarks[:, :, 1])
    plt.colorbar()
    plt.title('Camera 2 Dark')


med_ch1 = np.median(finalDarks[:, :, 0])
mean_ch1 = np.mean(finalDarks[:, :, 0])
med_ch2 = np.median(finalDarks[:, :, 1])
mean_ch2 = np.mean(finalDarks[:, :, 1])

print(' ')
print('Channel 1: median = %f, mean = %f' % (med_ch1, mean_ch1))
print('Channel 2: median = %f, mean = %f' % (med_ch2, mean_ch2))
print(' ')

if saveData:
    saveFilename = saveFilePref + filePref
    np.savez(saveFilename, allSummedDarks=allSummedIms, finalDarks=finalDarks, sigmaDarks=sigmaDarks)