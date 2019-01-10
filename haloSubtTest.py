import numpy as np
from scipy import ndimage
from scipy import io
from astropy.io import fits
import matplotlib.pyplot as plt
from copy import copy
import polzimtools as plz
import time

from scipy.optimize import curve_fit

# First get radial profile using CALIBRATOR star if star resolved:
# Get raw Stokes Images:

profLim = 100. #120 # Max no. of pixels to include in radial profile
radCut = 10. # Downweight region within this (i.e. the central peak).  ~10.
radCutSig = 1e4 # Inner region has its sigma set to this (elsewhere it's 1)

Ical=p.plotSingleIm(imType='I', cmap='viridis', norm=False)
Qcal=p.plotSingleIm(imType='Q', cmap='viridis', norm=False)
Ucal=p.plotSingleIm(imType='U', cmap='viridis', norm=False)


def radial_profile(data, center):
    x, y = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile

def gaus(x,a,x0,sigma, c=0):
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + c

def halofunc(x, a, b, c):
    return a*np.exp(-x**b/c**2)

img = Ical
center = np.unravel_index(img.argmax(), img.shape)
print center
#center = (-fitsFile[0].header['LBOUND2']+1, -fitsFile[0].header['LBOUND1']+1)
radProf = radial_profile(img, center)
#radProf = radProf / np.max(radProf)

radProf = radProf[0:profLim]
px = np.arange(profLim, dtype='float')

plt.clf()
plt.plot(radProf)
plt.yscale('log')
plt.xscale('log')
# plt.plot(px, halofunc(px, 3e2, 1, 5.5))

radCut=10
sigs = np.ones(len(px))
sigs[0:radCut] = radCutSig #Ignore central peak
popt, pcov = curve_fit(halofunc, px, radProf, p0=(1e6, 1., 1.), sigma=sigs)
plt.plot(px, halofunc(px, popt[0], popt[1], popt[2]))

print('Fitted params: a = %e, b = %f, c = %f' % (popt[0], popt[1], popt[2]))


# Make radial profile:
dim = len(img)
pxx = np.arange(dim) - center[1]
pxy = np.arange(dim) - center[0]
xv, yv = np.meshgrid(pxx, pxy)
rpix = np.sqrt(xv ** 2 + yv ** 2)

haloImage = halofunc(rpix, popt[0], popt[1], popt[2])
plt.imshow(np.log10(haloImage), interpolation='nearest', cmap='viridis')

