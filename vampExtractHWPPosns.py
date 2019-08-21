
dataPath = '/Volumes/BigVampData/201903/20190321/'
outfilename = 'vampHWPLog_20190321.txt'


filePref = 'ABAUr_01_20190226_750-50_EmptySlot_0'
fileSuf = '_cam1.fits'
hwpWaitTime = 5
hwpWaitTimeSafety = 1
nFiles = 32
nFiles = 2

# outfilename = './testOut.txt'

doAllInDir = True
excludeStr = 'linked'


import astropy.io.fits as fits
from datetime import datetime
from datetime import timedelta
import os


# Instead, do all files in a directory:
if doAllInDir:
    allFilenames = []
    allVampPrefs = []
    for r, d, f in os.walk(dataPath):
        for file in f:
            if ('.fits' in file) and ('cam1' in file) and (excludeStr not in r):
                allFilenames.append(os.path.join(r,file))
                allVampPrefs.append(file.split('_')[0])
    nFiles = len(allFilenames)



outfile = open(outfilename, 'w')
headerStr = 'HWP_ANGLE MOVE_START MOVE_END VAMP_PREFIX'
outfile.write(headerStr + '\n')

for f in range(nFiles):
    if doAllInDir:
        filename = allFilenames[f]
        vampPref = allVampPrefs[f]
    else:
        filename = dataPath + filePref + '%d' % f + fileSuf
        vampPref = filePref.split('_')[0]
    print(filename)
    try:
        hdr = fits.getheader(filename, 1)
    except:
        print("File ignored, couldn't load FITS extension header")
        continue
    hwpAngle = hdr['A0HWPVAL']
    startTime = hdr['UTSTTIME']
    endTime = hdr['UTNDTIME']
    # Assume the HWP starts moving hwpWaitTime+hwpWaitTimeSafety seconds before startTime.
    d = datetime.strptime(startTime, '%Y%m%dT%H%M%S')
    dt = timedelta(seconds = hwpWaitTime + hwpWaitTimeSafety)
    hwpMoveStart = d-dt
    hwpMoveStartStr = hwpMoveStart.strftime('%Y%m%dT%H%M%S')

    print('start time: ' + startTime)
    print('HWP move start time: '+ hwpMoveStartStr)
    print('HWP moved to angle: ' + str(hwpAngle))

    outStr = ('%05.2f ' % hwpAngle + hwpMoveStartStr + ' ' + startTime + ' ' + vampPref + '\n')
    print(outStr)
    print(' ')
    outfile.write(outStr)


outfile.close()