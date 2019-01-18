import numba
from scipy import ndimage
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg.linalg import multi_dot as md

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

def QUdiag(Q, U):
    Q = np.float(Q)
    U = np.float(U)
    A = np.array(( (1, 0, 0, 0),
                   (0, Q, 0, 0),
                   (0, 0, U, 0),
                   (0, 0, 0, 1)) )
    return A

@numba.jit(nopython=True)
def rotImPolz(inIm, theta):
    imsz = inIm.shape[0]
    ang = np.float(theta) / 180 * np.pi * 2
    A = np.array(( (1, 0, 0, 0),
                   (0, np.cos(ang), np.sin(ang), 0),
                   (0, -np.sin(ang), np.cos(ang), 0),
                   (0, 0, 0, 1)) )

    imOut = np.zeros((imsz, imsz ,4))
    for x in range(imsz):
        for y in range(imsz):
            curInVec = inIm[x,y,:]
            imOut[x, y, :] = np.dot(A, curInVec)
    return imOut


@numba.jit(nopython=True)
def stokesToPol(inIm):
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

    # This code is all written out longhand to make JIT straightforwards
    theta = 0.
    ang = -theta / 180 * np.pi * 2
    A1 = np.array(((1, 0, 0, 0),
                   (0, np.cos(ang), np.sin(ang), 0),
                   (0, -np.sin(ang), np.cos(ang), 0),
                   (0, 0, 0, 1)))
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
    A_0 = np.dot(A1, np.dot(A2, A3))

    theta = 90.
    ang = -theta / 180 * np.pi * 2
    A1 = np.array(((1, 0, 0, 0),
                   (0, np.cos(ang), np.sin(ang), 0),
                   (0, -np.sin(ang), np.cos(ang), 0),
                   (0, 0, 0, 1)))
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
    A_90 = np.dot(A1, np.dot(A2, A3))

    theta = 45.
    ang = -theta / 180 * np.pi * 2
    A1 = np.array(((1, 0, 0, 0),
                   (0, np.cos(ang), np.sin(ang), 0),
                   (0, -np.sin(ang), np.cos(ang), 0),
                   (0, 0, 0, 1)))
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
    A_45 = np.dot(A1, np.dot(A2, A3))

    theta = 135.
    ang = -theta / 180 * np.pi * 2
    A1 = np.array(((1, 0, 0, 0),
                   (0, np.cos(ang), np.sin(ang), 0),
                   (0, -np.sin(ang), np.cos(ang), 0),
                   (0, 0, 0, 1)))
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
    A_135 = np.dot(A1, np.dot(A2, A3))

    for x in range(imsz):
        for y in range(imsz):
            curInVec = inIm[x,y,:]
            imH[x, y] = np.dot(A_0, curInVec)[0]
            imV[x, y] = np.dot(A_90, curInVec)[0]
            imH45[x, y] = np.dot(A_45, curInVec)[0]
            imV45[x, y] = np.dot(A_135, curInVec)[0]

    allIms = (imH, imV, imH45, imV45)
    return allIms


@numba.jit(nopython=True)
def matIm(inIm, mat): # General function to multiply Stokes image by Mueller matrix
    imsz = inIm.shape[0]
    imOut = np.zeros((imsz, imsz, 4))
    for x in range(imsz):
        for y in range(imsz):
            curInVec = inIm[x, y, :]
            imOut[x, y, :] = np.dot(mat, curInVec)
    return imOut