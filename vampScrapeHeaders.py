# Scrape vampires headers from directory of files into a pandas dataframe

dataPath = '/Volumes/BigVampData/201903/20190321/'
outfile = './hdrDataOut_test1.csv'

dataPath = '/Volumes/BigVampData/'
outfile = './hdrDataOut_allOnBVD_02'


dataPath = '/Volumes/BigVampData/201902/20190226/'
outfile = './hdrDataOut_20190226all'
dataPath = '/Volumes/BigVampData/201810/20181022/'
outfile = './hdrDataOut_20181022all'



excludes = ['dark', 'flat', 'DARK', 'FLAT', 'HaDifferential', 'Analysis',
            'polzCal']



import csv
import numpy as np
import matplotlib.pyplot as plt
import polzimtools as plz
import pandas as pd
import os
import astropy.io.fits as fits


# Get a list of filenames to look at
allFilenames = []
allVampPrefs = []
for r, d, f in os.walk(dataPath):
    for file in f:
        print('Walking ' + r)
        fullFname = os.path.join(r, file)
        if ('.fits' in file) and \
                (not any(e in fullFname for e in excludes)):
            allFilenames.append(os.path.join(r, file))
            allVampPrefs.append(file.split('_')[0])

nFiles = len(allFilenames)
print('Total files found: %d' % nFiles)

data = pd.DataFrame(columns=['UTSTTIME', 'UTNDTIME', 'RA', 'DEC',
                             'PAP', 'PAD', 'IMRA', 'AIRMASS'])
allRows = []
displayCount = 0
displayInt = 10
# Go through files and extract headers
for filename in allFilenames:

    if displayCount == displayInt:
        print('Extracting header from '+filename)
        displayCount = 0
    else:
        displayCount = displayCount + 1

    try:
        hdr = fits.getheader(filename, 1)
    except:
        print("File ignored, couldn't load FITS header")
        continue

    try:
        UTSTTIME = hdr['UTSTTIME']
        UTNDTIME = hdr['UTNDTIME']
        RA = hdr['RA']
        DEC = hdr['DEC']
        PAP = hdr['PAP']
        PAD = hdr['PAD']
        IMRA = hdr['IMRA']
        AIRMASS = hdr['AIRMASS']
    except:
        print("File ignored, couldn't find all keywords.")
        continue

    curRow = {'UTSTTIME': UTSTTIME, 'UTNDTIME': UTNDTIME, 'RA': RA, 'DEC': DEC,
                             'PAP':PAP, 'PAD': PAD, 'IMRA': IMRA,
                            'AIRMASS': AIRMASS, 'FILENAME': filename}
    allRows.append(curRow)

hdrData = pd.DataFrame(allRows)

hdrData.to_csv(outfile+'.csv')
hdrData.to_hdf(outfile+'.h5', 'hdrData')
# To read:
# hdrData = pd.read_hdf(infile, 'hdrData')




