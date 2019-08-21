

infile = 'hdrDataOut_allOnBVD_02.h5'
stride=100
outfile = None

infile = 'hdrDataOut_20181022all.h5'
stride=1
outfile = 'hdrDataOut_20181022all_augmented'


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import astropy.units as u
import astropy.coordinates as ac
from astropy.time import Time as aTime
import time
from astroplan import Observer

hdrData = pd.read_hdf(infile, 'hdrData')
nTot = hdrData.shape[0]
rowsToUse = np.arange(0, nTot, stride)
selData = hdrData.loc[rowsToUse, :]


##### Remove problem rows
selData = selData[(selData['PAD'] < 360) & (selData['PAD'] > -360)]

# This will remove the IMR minus-sign-problem data...
selData = selData[selData['IMRA'] > 0]


##### Add new columns

# # from astropy.coordinates import SkyCoord  # High-level coordinates
# # from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
# # from astropy.coordinates import Angle, Latitude, Longitude  # Angles
# import astropy.units as u
# # from astropy.time import Time
# # from astropy.coordinates import SkyCoord, EarthLocation, AltAz
# import astropy.coordinates as ac
# from astropy.time import Time as aTime
# coords = ["1:12:43.2 +1:12:43"]
# c = ac.SkyCoord(coords, frame=FK4, unit=(u.deg, u.hourangle), obstime="J1992.21")
# Or
# c = ac.SkyCoord.from_name('M33')


# Get alt / az and airmass
obsLoc = ac.EarthLocation.of_site('Subaru Telescope')
subaruObserver = Observer(location=obsLoc, name="Subaru", timezone="US/Hawaii")
counter = 0
counterInt = 10
startTime = time.time()


for index, curRow in selData.iterrows():
    # curRow = selData.loc[0]
    # print(curRow)
    timeString = curRow['UTSTTIME']
    obsDT = datetime.strptime(timeString, '%Y%m%dT%H%M%S')
    obsTime = aTime(obsDT)

    raIn = curRow['RA']
    decIn = curRow['DEC']
    ra = raIn[0:2] + ' ' + raIn[3:5] + ' ' + raIn[6:]
    dec = decIn[0:3] + ' ' + decIn[4:6] + ' ' + decIn[7:]
    curStar = ac.SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

    altaz = curStar.transform_to(ac.AltAz(location=obsLoc, obstime=obsTime))
    airmass = altaz.secz

    selData.loc[curRow.name, 'ALT'] = altaz.alt.value
    selData.loc[curRow.name, 'AZ'] = altaz.az.value
    selData.loc[curRow.name, 'CALC_AIRMASS'] = airmass.value

    parang = subaruObserver.parallactic_angle(obsTime, curStar).value / np.pi * 180
    selData.loc[curRow.name, 'PARANG'] = parang

    if counter == counterInt:
        counter = 0
        print('%s iterations took %f seconds' % (counterInt, time.time()-startTime))
        startTime = time.time()
    else:
        counter = counter + 1

plt.plot(selData['ALT'], selData['IMRA'], '.')
plt.xlabel('ALT')
plt.ylabel('IMRA')

plt.plot(selData['ALT'], selData['PAD'], '.')
plt.xlabel('ALT')
plt.ylabel('PAD')

# plt.plot(selData['PAD'], selData['IMRA'], '.')

plt.plot(selData['PARANG'], selData['PAD'], '.')
plt.xlabel('PARANG')
plt.ylabel('PAD')


if outfile is not None:
    selData.to_csv(outfile + '.csv')
    selData.to_hdf(outfile + '.h5', 'hdrData')