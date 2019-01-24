import csv
import numpy as np
import matplotlib.pyplot as plt
import polzimtools as plz
import pandas as pd
from lookatCalTools import *

dataPath = '../SimplePDI_DATA/'
inFile = 'Measurements 750nm-Condensed A 20100124a.csv'


pd.set_option('display.max_columns', 100)

# with open(dataPath+inFile, 'r') as f:
#     reader = csv.reader(f)
#     inList = list(reader)

data = pd.read_csv(dataPath+inFile)
data = data.drop(columns='ID')

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


# Subset of cal0 data with scalar fractional polz:
# subs1 = data.loc[:, ['Target', 'IMR', 'I', 'Q', 'U', 'p_s', 'pQ_s', 'pU_s', 'True p', 'True pQ',
#                        'True pU', 'diff p', 'diff pQ', 'diff pU']]

subs1unpol = data[data['True p'] <= 1e-3]
subs1unpol = subs1unpol.sort_values('IMR')
subs1pol = data[data['True p'] > 1e-3]
subs1pol = subs1pol.sort_values('IMR')


colsToPlot = ['diff p', 'diff pQ', 'diff pU']
colsToPlot = ['True p', 'p_s', 'True pQ', 'pQ_s', 'True pU', 'pU_s']
# colsToPlot = ['True pQ', 'pQ_s', 'True pU', 'pU_s']


plotPolvals(subs1unpol, colsToPlot, fignum=1, ylims=(-0.06, 0.06))
plotPolvals(subs1pol, colsToPlot, fignum=2, ylims=(-0.12, 0.12))


colsToPlot = ['diff p 3a', 'diff pQ 3a', 'diff pU 3a']
colsToPlot = ['True p', 'p_s 3a', 'True pQ', 'pQ_s 3a', 'True pU', 'pU_s 3a']
# colsToPlot = ['True pQ', 'pQ_s 3a', 'True pU', 'pU_s 3a']

plotPolvals(subs1unpol, colsToPlot, fignum=3, ylims=(-0.06, 0.06))
plotPolvals(subs1pol, colsToPlot, fignum=4, ylims=(-0.12, 0.12))


colsToPlot = ['diff p 3b', 'diff pQ 3b', 'diff pU 3b']
colsToPlot = ['True p', 'p_s 3b', 'True pQ', 'pQ_s 3b', 'True pU', 'pU_s 3b']
# colsToPlot = ['True pQ', 'pQ_s 3b', 'True pU', 'pU_s 3b']

plotPolvals(subs1unpol, colsToPlot, fignum=5, ylims=(-0.06, 0.06))
plotPolvals(subs1pol, colsToPlot, fignum=6, ylims=(-0.12, 0.12))


#
# # Only include unpolarised stars
# inds = np.arange(len(subs1unpol))
# targets = subs1unpol.loc[:, 'Target'].values
# imrs = subs1unpol.loc[:, 'IMR'].values
# labels=[]
# for k in range(len(imrs)):
#     labels.append(targets[k] + ' (' + str(imrs[k]) + ')')
# plt.figure(5)
# plt.clf()
# ax = plt.gca()
# subs1unpol.loc[:, colsToPlot].plot.bar(ax=ax)
# plt.xticks(inds, labels, rotation=90)
# plt.ylim([-0.06, 0.06])
# plt.tight_layout()
#
# # Only include polarised stars
# inds = np.arange(len(subs1pol))
# targets = subs1pol.loc[:, 'Target'].values
# imrs = subs1pol.loc[:, 'IMR'].values
# labels=[]
# for k in range(len(imrs)):
#     labels.append(targets[k] + ' (' + str(imrs[k]) + ')')
# plt.figure(6)
# plt.clf()
# ax = plt.gca()
# subs1pol.loc[:, colsToPlot].plot.bar(ax=ax)
# plt.xticks(inds, labels, rotation=90)
# plt.ylim([-0.12, 0.12])
# plt.tight_layout()
#
#
#
