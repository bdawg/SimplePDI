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

#dataPath = '/Users/bnorris/DontBackup/preprocTestdata/'
#filePrefs = ['RYTau_01_20180106_750-50_EmptySlot_0',
#             'RYTau_01__RESTARTED___20180106_750-50_EmptySlot_0']
#nSubFilesAll = [4,4]
#outputPath = '/Users/bnorris/DontBackup/preprocTestdataOut/'
#outFilePref = 'RYTau_01_20180106_750_Fullpupil_'

# dataPath = '/import/pendragon1/snert/VAMPIRES/VAMPIRESData_201801/20180106/'
# dataPath = '/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201801/20180106/'
# filePrefs = ['HD283691_02_20180106_750-50_EmptySlot_0']
# nSubFilesAll = [16]
# outputPath = '/import/pendragon1/snert/VAMPIRES/VAMPIRESData_201801/compileData_Takayuki/HD283691_02/'
# outputPath = '/Volumes/PENDRAGON1/snert/VAMPIRES/VAMPIRESData_201801/compileData_Takayuki/HD283691_02/'
# outFilePref = 'HD283691_02_20180106_750_Fullpupil_'


def preprocVampiresFunc(dataPath, filePrefs, nSubFilesAll, outputPath, outFilePref):

    fileExtn = '.fits'
    startFileNum = 0

    # A list of keywords (for two extns) whose FITS header entries will be copied to the output file.
    copyKeywords0 = ['HIERARCH acqinttime', 'ACQNCYCS', 'HIERARCH acqflcoffset', 'HIERARCH acqstartdelay',
                     'HIERARCH acqstartdate', 'HIERARCH acqstarttime', 'HIERARCH acqdeltatime',
                     'HWP', 'DATE']
    copyKeywords1 = ['UTSTTIME', 'UTNDTIME', 'LOCTIME', 'LOOPN', 'TOTPOS', 'POSN', 'RA', 'DEC',
                     'EMGAIN', 'QWP1', 'QWP2', 'FILTER', 'MASK', 'MANGLE', 'AO188HWP', 'A0HWPVAL', 'PAP',
                     'PAD', 'IMRA', 'AIRMASS']

    ##############################################################################################
    import numpy as np
    from scipy import ndimage
    from scipy import io
    from astropy.io import fits
    from copy import copy
    import time

    # Read a FITS file to get sizes
    curFilenumStr = '%d' % startFileNum
    curCamStr = '_cam1'
    curFilename = dataPath + filePrefs[0] + curFilenumStr + curCamStr + fileExtn
    hdulist = fits.open(curFilename)
    curHDU = hdulist[0]
    curCube = np.transpose(curHDU.data)
    nFrms = curCube.shape[2]
    dim = curCube.shape[0]

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
                    curCamStrOut = '_1'
                else:
                    curCamStr = '_cam2'
                    curCamStrOut = '_2'

                curFilenumStr = '%d' % curFileNum
                curFilename = dataPath + filePref + curFilenumStr + curCamStr + fileExtn

                print('Reading file %s' % curFilename)
                hdulist = fits.open(curFilename)
                curHDU = hdulist[0]
                curSuperCube = np.transpose(curHDU.data)  # 'Super' since both FLC states

                # TODO - Get both headers (from both extns) and add to output file
                # TODO - also include original filename
                # # Get PA from header (use PAD, and assume instrument offset added later)
                # pa = hdulist[1].header['PAD']
                # allRawPAs.append(pa)

                # Now sort frames into the 2 FLC states and discard 1st 2 frames
                nSubFrms = nFrms // 2 - 1
                curCube_FLC1 = np.zeros((dim, dim, nSubFrms), dtype='uint16')
                curCube_FLC2 = np.zeros((dim, dim, nSubFrms), dtype='uint16')
                curFLC = 1
                count = 0
                for k in range(2, nFrms):
                    if curFLC == 1:
                        curCube_FLC1[:, :, count] = curSuperCube[:, :, k]
                    else:
                        curCube_FLC2[:, :, count] = curSuperCube[:, :, k]
                        count = count + 1

                    if curFLC == 1:
                        curFLC = 2
                    else:
                        curFLC = 1

                # Save output files
                for l in range(0, 2):
                    if l == 0:
                        curCube = curCube_FLC1
                        curFLCStr = '_A'
                    else:
                        curCube = curCube_FLC2
                        curFLCStr = '_B'

                    curOutFilenumStr = '%03d' % curOutFnum
                    curOutFilename = outputPath + outFilePref + curOutFilenumStr + curCamStrOut \
                            + curFLCStr + fileExtn
                    print('Writing file ' + curOutFilename)
                    #fits.writeto(curOutFilename, np.transpose(curCube))
                    data = np.transpose(curCube)
                    outHDU = fits.PrimaryHDU(data)
                    inHdr0 = hdulist[0].header
                    inHdr1 = hdulist[1].header
                    for keyword in copyKeywords0:
                        outHDU.header.append(inHdr0.cards[keyword])
                    for keyword in copyKeywords1:
                        outHDU.header.append(inHdr1.cards[keyword])
                    outHDU.writeto(curOutFilename)



            curOutFnum = curOutFnum + 1

























