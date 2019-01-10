import numpy as np


dataPath = "/Users/bnorris/DontBackup/simplePDIdata/Ha/"

filenames = ('allSummedImsCube_omiCet_01_HaDifferential-.npz',
             'allSummedImsCube_omiCet_02_HaDifferential-.npz',
             'allSummedImsCube_omiCet_04_HaDifferential-.npz')

joinedIms = np.zeros((256, 256, 0, 2, 2))
joinedPAs = np.zeros(0)

for filename in filenames:
    npzfile = np.load(dataPath+filename)
    allSummedIms = npzfile['arr_0']
    allPAs = npzfile['arr_1']
    joinedIms = np.concatenate((joinedIms, allSummedIms), axis=2)
    joinedPAs = np.concatenate((joinedPAs, allPAs))

saveFilename = dataPath + 'CombinedCube'
np.savez(saveFilename, joinedIms, joinedPAs)