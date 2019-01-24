import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



def plotPolvals(df, colsToPlot, fignum=1, ylims=(-0.1, 0.1)):
    inds = np.arange(len(df))
    targets = df.loc[:, 'Target'].values
    imrs = df.loc[:, 'IMR'].values
    labels = []
    for k in range(len(imrs)):
        labels.append(targets[k] + ' (' + str(imrs[k]) + ')')
    plt.figure(fignum)
    plt.clf()
    ax = plt.gca()
    df.loc[:, colsToPlot].plot.bar(ax=ax)
    plt.xticks(inds, labels, rotation=90)
    plt.ylim(ylims)
    plt.tight_layout()