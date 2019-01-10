from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg.linalg import multi_dot as md
from time import clock as cl
from scipy import ndimage

infile = 'toyDisk.fits'



# Useful Mueller matrices:
def linPol(gamma=0):
    # gamma is in degrees
    ang = np.float(gamma)/180*np.pi * 2
    A = np.array([ [1, np.cos(ang), 0, 0],
                    [np.cos(ang), 1, 0, 0],
                    [0, 0, np.sin(ang), 0],
                    [0, 0, 0, np.sin(ang)] ])
    return A

def rot(theta):
    # theta is in degrees
    ang = np.float(theta)/180*np.pi * 2
    A = np.array([[1, 0, 0, 0],
                  [0, np.cos(ang), np.sin(ang), 0],
                  [0, -np.sin(ang), np.cos(ang), 0],
                  [0, 0, 0, 1] ])
    return A

def rotLP(theta, gamma=0):
    #t = time.clock()
    A = md([rot(-theta), linPol(gamma), rot(theta)])
    #print time.clock() - t
    #A = np.zeros([4,4])
    return A




t = cl()
hdulist = fits.open(infile)
curHDU = hdulist[0]
inIm = np.transpose(curHDU.data, (1, 2, 0))
imsz = inIm.shape[0]
imI = inIm[:,:,0]
imQ = inIm[:,:,1]
imU = inIm[:,:,2]
imV = inIm[:,:,3]

imH = np.zeros([imsz,imsz])
imV = np.zeros([imsz,imsz])
imH45 = np.zeros([imsz,imsz])
imV45 = np.zeros([imsz,imsz])

for x in range(imsz):
    for y in range(imsz):
        curInVec = inIm[x,y,:]
        t = cl()
        theta = 0.
        imH[x,y] = np.dot(rotLP(theta), curInVec)[0]
        theta = 90.
        imV[x,y] = np.dot(rotLP(theta), curInVec)[0]
        theta = 45.
        imH45[x,y] = np.dot(rotLP(theta), curInVec)[0]
        theta = 135.
        imV45[x,y] = np.dot(rotLP(theta), curInVec)[0]
        print cl()-t
allIms = [imH, imV, imH45, imV45]
#return allIms



t = cl()
##allIms = makeIms(infile)
print cl()-t


print 'Done.'


