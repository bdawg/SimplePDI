import csv
import numpy as np
import matplotlib.pyplot as plt
import polzimtools as plz
import pandas as pd
from lookatCalTools import *

desired_width = 320
pd.set_option('display.width', desired_width)

dataPath = '../SimplePDI_DATA/'
# inFile = 'Measurements 750nm-Condensed A 20100124a.csv'
# inFile = 'Measurements 750nm REVISION 2-Condensed A 20100124a.csv'
# inFile = 'Measurements 750nm-Condensed A aper8.csv'
inFile = 'Measurements 750nm REVISION 2-Condensed A_correctFactorof2.csv'
inFile = 'Measurements 750nm REVISION 2-Condensed A withMaxpix.csv'
# inFile = 'Measurements 750nm REV2 NO POLZROT-Condensed A_correctFactorof2.csv'

# inFile = 'Meas 750nm R2 MinusPA WIth+54-Condensed A.csv'
# inFile = 'Meas 750nm R2 MinusPA-Condensed A.csv'

inFile = 'M750R2 - Pol+Zero cals grouped-Condensed A.csv'



sortByIMR = True

figsize1 = (7, 5)
figsize2 = (4, 3)

pd.set_option('display.max_columns', 100)

# with open(dataPath+inFile, 'r') as f:
#     reader = csv.reader(f)
#     inList = list(reader)

data = pd.read_csv(dataPath+inFile)
# data = data.drop(columns='ID')


# Drop rows where data is saturated
satThresh = 30000
# data = data[data['Abs max pix'] <= satThresh]




# Make new columns showing difference from true values, using scalar measurements
# data['True p'] = np.sqrt(data['True pQ']**2 +data['True pU']**2)
data['diff p'] = data['p_s'] - data['True p']
data['diff pQ'] = data['pQ_s'] - data['True pQ']
data['diff pU'] = data['pU_s'] - data['True pU']
data['diff p 3a'] = data['p_s 3a'] - data['True p']
data['diff pQ 3a'] = data['pQ_s 3a'] - data['True pQ']
data['diff pU 3a'] = data['pU_s 3a'] - data['True pU']
data['diff p 3b'] = data['p_s 3b'] - data['True p']
data['diff pQ 3b'] = data['pQ_s 3b'] - data['True pQ']
data['diff pU 3b'] = data['pU_s 3b'] - data['True pU']



# Make new measurements rotated by specified offset
offsetAngle = 54.
# offsetAngle = 0
pQ_newVals = []
pU_newVals = []
for index, row in data.iterrows():
    pQ = row['pQ_s']
    pU = row['pU_s']
    theta1 = 0.5 * np.arctan(pU / pQ) / np.pi * 180
    theta2 = 0.5 * np.arctan2(pU, pQ) / np.pi * 180
    p_calc = np.sqrt(pQ ** 2 + pU ** 2)

    # print('Theta 1: %f' % theta1)
    # print('Theta 2: %f' % theta2)
    # print('Calcd p: %f' % p_calc)

    newTheta = theta2 + offsetAngle
    pQ_new = p_calc * np.cos(2 * newTheta / 180 * np.pi)
    pU_new = p_calc * np.sin(2 * newTheta / 180 * np.pi)
    # print('new pQ: %f' % pQ_new)
    # print('new pU: %f' % pU_new)
    # print('new theta: %f' % newTheta)
    # print(' ')
    pQ_newVals.append(pQ_new)
    pU_newVals.append(pU_new)

data['pQ_s'] = pQ_newVals
data['pU_s'] = pU_newVals
data['pQ_s Rotated'] = pQ_newVals
data['pU_s Rotated'] = pU_newVals



# Subset of cal0 data with scalar fractional polz:
# subs1 = data.loc[:, ['Target', 'IMR', 'I', 'Q', 'U', 'p_s', 'pQ_s', 'pU_s', 'True p', 'True pQ',
#                        'True pU', 'diff p', 'diff pQ', 'diff pU']]

subs1unpol = data[data['True p'] <= 1e-3]
subs1pol = data[data['True p'] > 1e-3]
if sortByIMR:
    subs1unpol = subs1unpol.sort_values('IMR')
    subs1pol = subs1pol.sort_values('IMR')




colsToPlot = ['diff p', 'diff pQ', 'diff pU']
colsToPlot = ['True p', 'p_s', 'True pQ', 'pQ_s', 'True pU', 'pU_s']
colsToPlot = ['True pQ', 'pQ_s', 'True pU', 'pU_s']

plt.figure(1, figsize=figsize1)
plotPolvals(subs1unpol, colsToPlot, fignum=1, ylims=(-0.06, 0.06))
plt.figure(2, figsize=figsize1)
plotPolvals(subs1pol, colsToPlot, fignum=2, ylims=(-0.12, 0.12))
# colsToPlot = ['diff p 3a', 'diff pQ 3a', 'diff pU 3a']
# colsToPlot = ['True p', 'p_s 3a', 'True pQ', 'pQ_s 3a', 'True pU', 'pU_s 3a']
# colsToPlot = ['True pQ', 'pQ_s 3a', 'True pU', 'pU_s 3a']
#
# colsToPlot = ['True pQ', 'pQ_s Rotated', 'True pU', 'pU_s Rotated']
#
# plotPolvals(subs1unpol, colsToPlot, fignum=3, ylims=(-0.06, 0.06))
# plotPolvals(subs1pol, colsToPlot, fignum=4, ylims=(-0.12, 0.12))


# colsToPlot = ['diff p 3b', 'diff pQ 3b', 'diff pU 3b']
# colsToPlot = ['True p', 'p_s 3b', 'True pQ', 'pQ_s 3b', 'True pU', 'pU_s 3b']
# colsToPlot = ['True pQ', 'pQ_s 3b', 'True pU', 'pU_s 3b']
#
# plotPolvals(subs1unpol, colsToPlot, fignum=5, ylims=(-0.06, 0.06))
# plotPolvals(subs1pol, colsToPlot, fignum=6, ylims=(-0.12, 0.12))


# Plot Stokes vector for unpolarised sources
lim = (-0.07, 0.07)
dataToPlot = subs1unpol
plt.figure(3, figsize=figsize2)
plt.clf()
n = len(dataToPlot.index)
colors = iter(plt.cm.rainbow(np.linspace(0,1,n)))
for index, row in dataToPlot.iterrows():
    pQ = row['pQ_s Rotated']
    pU = row['pU_s Rotated']
    c = next(colors)
    plt.plot((0, pQ), (0, pU), '-o', c=c, linewidth=1, markersize=4)
plt.xlim(lim)
plt.ylim(lim)


# Plot Stokes vector for polarised sources
lim = (-0.07, 0.07)
dataToPlot = subs1pol
colorDict = {'HD 43384': 'r', 'HD 154445': 'm', 'HD 155528': 'b', 'HD 198478': 'c', 'HD 183143': 'g'}
plt.figure(4, figsize=figsize2)
plt.clf()
n = len(dataToPlot.index)
colors = iter(plt.cm.rainbow(np.linspace(0,1,n)))
for index, row in dataToPlot.iterrows():
    pQ = row['pQ_s Rotated']
    pU = row['pU_s Rotated']
    targ = row['Target']
    formatStr = '-o'+colorDict[targ]
    plt.plot((0, pQ), (0, pU), formatStr, linewidth=1, markersize=4)

    pQ = row['True pQ']
    pU = row['True pU']
    formatStr = 'X' + colorDict[targ]
    plt.plot((pQ), (pU), formatStr, linewidth=1, markersize=8)
plt.xlim(lim)
plt.ylim(lim)

print(eval("subs1pol[['ID', 'Target', 'True pQ', 'True pU', 'pQ_s', 'pU_s', 'pQ_s Rotated', 'pU_s Rotated']]"))



##### Do calibration by differences
src2cal = [ [8, 6],
            [10, 2],
            [19, 2],
            [20, 1],
            [26, 7],
            [30, 7],
            [31, 7],
            [34, 48],
            [37, 50]
            ]

for pair in src2cal:
    srcRow = data[data['ID'] == pair[0]]
    calRow = data[data['ID'] == pair[1]]
    pQ = float(srcRow['pQ_s']) - float(calRow['pQ_s'])
    pU = float(srcRow['pU_s']) - float(calRow['pU_s'])
    data.loc[data['ID'] == pair[0], 'pQ_s cald'] = pQ
    data.loc[data['ID'] == pair[0], 'pU_s cald'] = pU
    data.loc[data['ID'] == pair[0], 'p_s cald'] = np.sqrt(pQ**2 + pU**2)

# dataToPlot = data[~np.isnan(data['pQ_s cald'])]

subs1pol = data[data['True p'] > 1e-3]
if sortByIMR:
    subs1pol = subs1pol.sort_values('IMR')

# Plot Stokes vector for polarised sources
lim = (-0.07, 0.07)
dataToPlot = subs1pol
colorDict = {'HD 43384': 'r', 'HD 154445': 'm', 'HD 155528': 'b', 'HD 198478': 'c', 'HD 183143': 'g'}
plt.figure(5, figsize=figsize2)
plt.clf()
n = len(dataToPlot.index)
colors = iter(plt.cm.rainbow(np.linspace(0,1,n)))
for index, row in dataToPlot.iterrows():
    pQ = row['pQ_s cald']
    pU = row['pU_s cald']
    targ = row['Target']
    formatStr = '-o'+colorDict[targ]
    plt.plot((0, pQ), (0, pU), formatStr, linewidth=1, markersize=4)

    pQ = row['True pQ']
    pU = row['True pU']
    formatStr = 'X' + colorDict[targ]
    plt.plot((pQ), (pU), formatStr, linewidth=1, markersize=8)
plt.xlim(lim)
plt.ylim(lim)

print(eval("subs1pol[['ID', 'Target', 'True p', 'True pQ', 'True pU', 'p_s', 'pQ_s', 'pU_s', "
           "'p_s cald', 'pQ_s cald', 'pU_s cald']]"))











