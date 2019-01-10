from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg.linalg import multi_dot as md
from time import clock as cl
import timeit
import numba
from scipy import ndimage

infile = 'toyDisk.fits'
plt.ioff()



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




def makeIms(infile):
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
            #t = cl()
            theta = 0.
            imH[x,y] = np.dot(rotLP(theta), curInVec)[0]
            theta = 90.
            imV[x,y] = np.dot(rotLP(theta), curInVec)[0]
            theta = 45.
            imH45[x,y] = np.dot(rotLP(theta), curInVec)[0]
            theta = 135.
            imV45[x,y] = np.dot(rotLP(theta), curInVec)[0]
            #print cl()-t
    allIms = [imH, imV, imH45, imV45]
    return allIms


#@numba.jit(nopython=True, locals={rot: numba.types.float64[:,:], linPol: numba.types.float64[:,:]})
@numba.jit(nopython=True)
def makeImsJit(inIm):
    imsz = inIm.shape[0]
    # imI = inIm[:,:,0]
    # imQ = inIm[:,:,1]
    # imU = inIm[:,:,2]
    # imV = inIm[:,:,3]

    imH = np.zeros((imsz,imsz))
    imV = np.zeros((imsz,imsz))
    imH45 = np.zeros((imsz,imsz))
    imV45 = np.zeros((imsz,imsz))
    gamma = 0.

    for x in range(imsz):
        for y in range(imsz):
            curInVec = inIm[x,y,:]
            #t = cl()

            theta = 0.
            ang = -theta / 180 * np.pi * 2
            A1 = np.array(( (1, 0, 0, 0),
                         (0, np.cos(ang), np.sin(ang), 0),
                         (0, -np.sin(ang), np.cos(ang), 0),
                         (0, 0, 0, 1) ))
            ang = gamma / 180 * np.pi * 2
            A2 = np.array(((1, np.cos(ang), 0, 0),
                          (np.cos(ang), 1, 0, 0),
                          (0, 0, np.sin(ang), 0),
                          (0, 0, 0, np.sin(ang))))
            ang = theta / 180 * np.pi * 2
            A3 = np.array(((1, 0, 0, 0),
                          (0, np.cos(ang), np.sin(ang), 0),
                          (0, -np.sin(ang), np.cos(ang), 0),
                          (0, 0, 0, 1)))
            A = np.dot(A1, np.dot(A2, A3))
            imH[x, y] = np.dot(A, curInVec)[0]
            # p = np.dot(A, curInVec)[0]
            # imH.itemset((x,y), p)

            theta = 90.
            ang = -theta / 180 * np.pi * 2
            A1 = np.array(( (1, 0, 0, 0),
                         (0, np.cos(ang), np.sin(ang), 0),
                         (0, -np.sin(ang), np.cos(ang), 0),
                         (0, 0, 0, 1) ))
            ang = gamma / 180 * np.pi * 2
            A2 = np.array(((1, np.cos(ang), 0, 0),
                          (np.cos(ang), 1, 0, 0),
                          (0, 0, np.sin(ang), 0),
                          (0, 0, 0, np.sin(ang))))
            ang = theta / 180 * np.pi * 2
            A3 = np.array(((1, 0, 0, 0),
                          (0, np.cos(ang), np.sin(ang), 0),
                          (0, -np.sin(ang), np.cos(ang), 0),
                          (0, 0, 0, 1)))
            A = np.dot(A1, np.dot(A2, A3))
            imV[x, y] = np.dot(A, curInVec)[0]
            # p = np.dot(A, curInVec)[0]
            # imV.itemset((x, y), p)

            theta = 45.
            ang = -theta / 180 * np.pi * 2
            A1 = np.array(( (1, 0, 0, 0),
                         (0, np.cos(ang), np.sin(ang), 0),
                         (0, -np.sin(ang), np.cos(ang), 0),
                         (0, 0, 0, 1) ))
            ang = gamma / 180 * np.pi * 2
            A2 = np.array(((1, np.cos(ang), 0, 0),
                          (np.cos(ang), 1, 0, 0),
                          (0, 0, np.sin(ang), 0),
                          (0, 0, 0, np.sin(ang))))
            ang = theta / 180 * np.pi * 2
            A3 = np.array(((1, 0, 0, 0),
                          (0, np.cos(ang), np.sin(ang), 0),
                          (0, -np.sin(ang), np.cos(ang), 0),
                          (0, 0, 0, 1)))
            A = np.dot(A1, np.dot(A2, A3))
            imH45[x, y] = np.dot(A, curInVec)[0]
            # p = np.dot(A, curInVec)[0]
            # imH45.itemset((x, y), p)

            theta = 135.
            ang = -theta / 180 * np.pi * 2
            A1 = np.array(( (1, 0, 0, 0),
                         (0, np.cos(ang), np.sin(ang), 0),
                         (0, -np.sin(ang), np.cos(ang), 0),
                         (0, 0, 0, 1) ))
            ang = gamma / 180 * np.pi * 2
            A2 = np.array(((1, np.cos(ang), 0, 0),
                          (np.cos(ang), 1, 0, 0),
                          (0, 0, np.sin(ang), 0),
                          (0, 0, 0, np.sin(ang))))
            ang = theta / 180 * np.pi * 2
            A3 = np.array(((1, 0, 0, 0),
                          (0, np.cos(ang), np.sin(ang), 0),
                          (0, -np.sin(ang), np.cos(ang), 0),
                          (0, 0, 0, 1)))
            A = np.dot(A1, np.dot(A2, A3))
            imV45[x, y] = np.dot(A, curInVec)[0]
            # p = np.dot(A, curInVec)[0]
            # imV45.itemset((x, y), p)

    allIms = (imH, imV, imH45, imV45)
    return allIms



hdulist = fits.open(infile)
curHDU = hdulist[0]
inIm = np.transpose(curHDU.data, (1, 2, 0))
inIm = inIm.astype('f8')

# t = cl()
# allIms = makeImsJit(inIm)
# print cl()-t


# Ok. Make a fake set of data, where the 'true' (Stokes) image is rotated, then
# 'analysed' with the linear polarisers.

paRange = range(-45, 45, 10)
allRotatedImages = []
imsz = inIm.shape[0]
for angle in paRange:
    print "Doing angle %f" % angle
    curStokesIm = np.zeros((imsz, imsz, 4))
    for s in range(4):
        curStokesIm[:,:,s] = ndimage.interpolation.rotate(inIm[:,:,s],angle,reshape=False)
    curPolIm = makeImsJit(curStokesIm)
    allRotatedImages.append(curPolIm)

    plt.imshow(curPolIm[0])
    plt.pause(0.001)

print 'Done.'


