# -*- coding: utf-8 -*-

# -----
# Setup
# -----

# script

# import statements
import sys
import os
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import QFileDialog, QApplication
import librosa
import dsp.filterFuncs
from scipy import stats
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns
import statsmodels.api as sm

# set plot parameters

sns.set_style('white')
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
mpl.rcParams.update({'font.size': 16})
mpl.rcParams['figure.autolayout'] = True
mpl.rcParams['mathtext.fontset'] = 'stix'

SMALL_SIZE = 9
MEDIUM_SIZE = 12
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE,
       labelsize=MEDIUM_SIZE)    # fontsize of the axes title and x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

mycolours = [(0, 102, 255), (0, 204, 153), (255, 0, 102), (74, 111, 152),
             (251, 164, 49), (204, 153, 255), (90, 192, 255), (90, 255, 243),
             (255, 90, 192), (164, 201, 242), (255, 254, 139), (255, 243, 255)]
mycolours = [tuple(shade/255 for shade in colour) for colour in mycolours]


# enable copy-on-write mode for Pandas (will be default from Pandas 3.0)
pd.options.mode.copy_on_write = True

# import data
# Part A
app = QApplication(sys.argv)
fileExts = "*.csv"
partADataFilePath = list(QFileDialog.getOpenFileName(filter=fileExts))[0]
partADataBySubj = pd.read_csv(partADataFilePath, index_col=False)

# Part B
partBDataFilePath = list(QFileDialog.getOpenFileName(filter=fileExts))[0]
partBDataBySubj = pd.read_csv(partBDataFilePath, index_col=False)

# data by stimulus
dataByStimFilePath = list(QFileDialog.getOpenFileName(filter=fileExts))[0]
dataByStim = pd.read_csv(dataByStimFilePath, index_col=0)
dataByStimAFilePath = list(QFileDialog.getOpenFileName(filter=fileExts))[0]
dataByStimTestA = pd.read_csv(dataByStimAFilePath, index_col=0)
dataByStimBFilePath = list(QFileDialog.getOpenFileName(filter=fileExts))[0]
dataByStimTestB = pd.read_csv(dataByStimBFilePath, index_col=0)

# pre and post test data
preTestDataFilePath = list(QFileDialog.getOpenFileName(filter=fileExts))[0]
preTestResponses = pd.read_csv(preTestDataFilePath, index_col=0)
postTestDataFilePath = list(QFileDialog.getOpenFileName(filter=fileExts))[0]
postTestResponses = pd.read_csv(postTestDataFilePath, index_col=0)

# categorise columns
partADataBySubj['SNRlevel'] = pd.Categorical(partADataBySubj['SNRlevel'],
                                             ["No UAS", "-16", "-10", "-4",
                                              "2", "8"])
partADataBySubj['UASLAeq'] = pd.Categorical(partADataBySubj['UASLAeq'],
                                            ["No UAS", "42", "48", "54", "60"])
partADataBySubj['UASOperation'] = pd.Categorical(partADataBySubj['UASOperation'],
                                                 ["No UAS", "Landing",
                                                  "Flyby", "Takeoff"])


partBDataBySubj['SNRlevel'] = pd.Categorical(partBDataBySubj['SNRlevel'],
                                             ["No UAS", "2", "8"])
partBDataBySubj['UASLAeq'] = pd.Categorical(partBDataBySubj['UASLAeq'],
                                            ["No UAS", "54", "60"])
partBDataBySubj['UASOperation'] = pd.Categorical(partBDataBySubj['UASOperation'],
                                                 ["No UAS", "Flyby"])

dataByStimTestA['SNRlevel'] = pd.Categorical(dataByStimTestA['SNRlevel'],
                                             ["No UAS", "-16", "-10", "-4",
                                              "2", "8"])
dataByStimTestA['UASLAeq'] = pd.Categorical(dataByStimTestA['UASLAeq'],
                                            ["No UAS", "42", "48", "54", "60"])
dataByStimTestA['UASOperation'] = pd.Categorical(dataByStimTestA['UASOperation'],
                                                 ["No UAS", "Landing",
                                                  "Flyby", "Takeoff"])

dataByStimTestB['SNRlevel'] = pd.Categorical(dataByStimTestB['SNRlevel'],
                                             ["No UAS", "2", "8"])
dataByStimTestB['UASLAeq'] = pd.Categorical(dataByStimTestB['UASLAeq'],
                                            ["No UAS", "54", "60"])
dataByStimTestB['UASOperation'] = pd.Categorical(dataByStimTestB['UASOperation'],
                                                 ["No UAS", "Flyby"])


# --------------------
# Exploratory analysis
# --------------------

# Part A
# ------

# UAS LAeq

UASLAeqCats = list(partADataBySubj['UASLAeq'].sort_values().unique())

# Annoyance

fig, ax = plt.subplots(figsize=(7, 4))

y_data = [partADataBySubj[partADataBySubj['UASLAeq']
                          == LAeq]['Annoyance'].values
          for LAeq in UASLAeqCats]

xjitter = 0.03
yjitter = 0.06
x_data = [np.array([i] * len(d)) for i, d in enumerate(y_data)]
x_jittered = [x + stats.t(df=6, scale=xjitter).rvs(len(x)) for x in x_data]
y_jittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y)) for y in y_data]

violins = ax.violinplot(y_data,
                        positions=range(0, len(UASLAeqCats)),
                        widths=0.45,
                        bw_method='scott',
                        showmeans=False,
                        showmedians=False,
                        showextrema=False)
for ii, pc in enumerate(violins["bodies"]):
    pc.set_facecolor(mycolours[ii])
    pc.set_edgecolor([0.25, 0.25, 0.25])
    pc.set_linewidth(1)
    pc.set_alpha(0.25)

medianprops = dict(linewidth=4,
                   color=[0.2, 0.2, 0.2],
                   solid_capstyle="butt")

boxprops = dict(linewidth=2,
                color=[0.2, 0.2, 0.2])

ax.boxplot(y_data,
           positions=range(0, len(UASLAeqCats)),
           showfliers=False,
           showcaps=False,
           medianprops=medianprops,
           whiskerprops=boxprops,
           boxprops=boxprops,
           widths=0.25)

# Add jittered dots ----------------------------------------------
for x, y, color in zip(x_jittered, y_jittered, mycolours[0:len(UASLAeqCats)]):
    ax.scatter(x, y, s=10, color=color, alpha=0.2)

ax.set(yticks=range(0, 11), xticks=range(0, len(UASLAeqCats)),
       xticklabels=UASLAeqCats,
       xlabel=r"UAS $L_\mathrm{Aeq}$, dB",
       ylabel="Annoyance rating", ylim=[-0.5, 10.5])
plt.show()

# segregated by ambient environment
fig, ax = plt.subplots(figsize=(7, 4.65))

data = partADataBySubj.loc[:, ['AmbientEnv', 'UASLAeq', 'Annoyance']]
data = data.sort_values('UASLAeq')

sns.violinplot(data=data, y='Annoyance', split=True, hue='AmbientEnv',
               x='UASLAeq', cut=0, palette=mycolours[0:2], inner='quart',
               width=0.5, bw_method='scott')

plt.setp(ax.collections, alpha=0.3)
# Add jittered dots ----------------------------------------------
y_subdataPark = list()
y_subdataStreet = list()
for ii in range(0, len(UASLAeqCats)):
    y_subdataPark.append(partADataBySubj.loc[np.logical_and(
                                             partADataBySubj['AmbientEnv']
                                             == 'Park',
                                             partADataBySubj['UASLAeq']
                                             == UASLAeqCats[ii]),
                                             'Annoyance'].values)
    y_subdataStreet.append(partADataBySubj.loc[np.logical_and(
                                               partADataBySubj['AmbientEnv']
                                               == 'Street',
                                               partADataBySubj['UASLAeq']
                                               == UASLAeqCats[ii]),
                                               'Annoyance'].values)
x_subdataPark = [np.array([i] * len(d)) for i, d in enumerate(y_subdataPark)]
x_subdataStreet = [np.array([i] * len(d))
                   for i, d in enumerate(y_subdataStreet)]
x_Parkjittered = [x - 0.05 + stats.t(df=6, scale=xjitter).rvs(len(x))
                  for x in x_subdataPark]
y_Parkjittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y))
                  for y in y_subdataPark]
x_Streetjittered = [x + 0.05 + stats.t(df=6, scale=xjitter).rvs(len(x))
                    for x in x_subdataStreet]
y_Streetjittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y))
                    for y in y_subdataStreet]
for x, y in zip(x_Parkjittered, y_Parkjittered):
    ax.scatter(x, y, s=10, color=mycolours[0], alpha=0.2)
for x, y in zip(x_Streetjittered, y_Streetjittered):
    ax.scatter(x, y, s=10, color=mycolours[1], alpha=0.2)

ax.set(yticks=range(0, 11), xticks=range(0, len(UASLAeqCats)),
       xticklabels=UASLAeqCats, xlabel=r"UAS $L_\mathrm{Aeq}$, dB",
       ylabel="Annoyance rating", ylim=[-0.5, 10.5])
ax.legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=2, fontsize=11)
plt.show()

# Valence & arousal
# segregated by ambient environment
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7, 7))

data0 = partADataBySubj.loc[:, ['AmbientEnv', 'UASLAeq', 'Valence']]
data0 = data0.sort_values('UASLAeq')
data1 = partADataBySubj.loc[:, ['AmbientEnv', 'UASLAeq', 'Arousal']]
data1 = data1.sort_values('UASLAeq')


sns.violinplot(ax=axs[0], data=data0, y='Valence', split=True, hue='AmbientEnv',
               x='UASLAeq', cut=0, palette=mycolours[0:2], inner='quart',
               width=0.5, bw_method=0.55, legend=True)
sns.violinplot(ax=axs[1], data=data1, y='Arousal', split=True, hue='AmbientEnv',
               x='UASLAeq', cut=0, palette=mycolours[0:2], inner='quart',
               width=0.5, bw_method=0.55, legend=False)

plt.setp(axs[0].collections, alpha=0.3)
plt.setp(axs[1].collections, alpha=0.3)
# Add jittered dots ----------------------------------------------
y_subdata0Park = list()
y_subdata0Street = list()
y_subdata1Park = list()
y_subdata1Street = list()
for ii in range(0, len(UASLAeqCats)):
    y_subdata0Park.append(partADataBySubj.loc[np.logical_and(
                                              partADataBySubj['AmbientEnv']
                                              == 'Park',
                                              partADataBySubj['UASLAeq']
                                              == UASLAeqCats[ii]),
                                              'Valence'].values)
    y_subdata0Street.append(partADataBySubj.loc[np.logical_and(
                                                partADataBySubj['AmbientEnv']
                                                == 'Street',
                                                partADataBySubj['UASLAeq']
                                                == UASLAeqCats[ii]),
                                                'Valence'].values)
    y_subdata1Park.append(partADataBySubj.loc[np.logical_and(
                                              partADataBySubj['AmbientEnv']
                                              == 'Park',
                                              partADataBySubj['UASLAeq']
                                              == UASLAeqCats[ii]),
                                              'Arousal'].values)
    y_subdata1Street.append(partADataBySubj.loc[np.logical_and(
                                                partADataBySubj['AmbientEnv']
                                                == 'Street',
                                                partADataBySubj['UASLAeq']
                                                == UASLAeqCats[ii]),
                                                'Arousal'].values)
x_subdata0Park = [np.array([i] * len(d)) for i, d in enumerate(y_subdata0Park)]
x_subdata0Street = [np.array([i] * len(d))
                    for i, d in enumerate(y_subdata0Street)]
x0_Parkjittered = [x - 0.05 + stats.t(df=6, scale=xjitter).rvs(len(x))
                   for x in x_subdata0Park]
y0_Parkjittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y))
                   for y in y_subdata0Park]
x0_Streetjittered = [x + 0.05 + stats.t(df=6, scale=xjitter).rvs(len(x))
                     for x in x_subdata0Street]
y0_Streetjittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y))
                     for y in y_subdata0Street]
for x, y in zip(x0_Parkjittered, y0_Parkjittered):
    axs[0].scatter(x, y, s=5, color=mycolours[0], alpha=0.2)
for x, y in zip(x0_Streetjittered, y0_Streetjittered):
    axs[0].scatter(x, y, s=5, color=mycolours[1], alpha=0.2)
x_subdata1Park = [np.array([i] * len(d)) for i, d in enumerate(y_subdata1Park)]
x_subdata1Street = [np.array([i] * len(d))
                    for i, d in enumerate(y_subdata1Street)]
x1_Parkjittered = [x - 0.05 + stats.t(df=6, scale=xjitter).rvs(len(x))
                   for x in x_subdata1Park]
y1_Parkjittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y))
                   for y in y_subdata1Park]
x1_Streetjittered = [x + 0.05 + stats.t(df=6, scale=xjitter).rvs(len(x))
                     for x in x_subdata1Street]
y1_Streetjittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y))
                     for y in y_subdata1Street]
for x, y in zip(x1_Parkjittered, y1_Parkjittered):
    axs[1].scatter(x, y, s=5, color=mycolours[0], alpha=0.2)
for x, y in zip(x1_Streetjittered, y1_Streetjittered):
    axs[1].scatter(x, y, s=5, color=mycolours[1], alpha=0.2)

axs[0].set(yticks=range(1, 6), xticks=range(0, len(UASLAeqCats)),
           xticklabels=UASLAeqCats, xlabel=None,
           ylabel="Valence rating", ylim=[0.5, 5.5])
axs[1].set(yticks=range(1, 6), xticks=range(0, len(UASLAeqCats)),
           xticklabels=UASLAeqCats, xlabel=r"UAS $L_\mathrm{Aeq}$, dB",
           ylabel="Arousal rating", ylim=[0.5, 5.5])
axs[0].legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=2,
              fontsize=11)
plt.show()


# UAS SNR

fig, ax = plt.subplots(figsize=(8.15, 4))

SNRCats = list(partADataBySubj['SNRlevel'].sort_values().unique())

y_data = [partADataBySubj[partADataBySubj['SNRlevel']
                          == SNR]['Annoyance'].values
          for SNR in SNRCats]

violins = ax.violinplot(y_data,
                        positions=range(0, len(SNRCats)),
                        widths=0.45,
                        bw_method='scott',
                        showmeans=False,
                        showmedians=False,
                        showextrema=False)
for ii, pc in enumerate(violins["bodies"]):
    pc.set_facecolor(mycolours[ii])
    pc.set_edgecolor([0.25, 0.25, 0.25])
    pc.set_linewidth(1)
    pc.set_alpha(0.25)

medianprops = dict(linewidth=4,
                   color=[0.2, 0.2, 0.2],
                   solid_capstyle="butt")

boxprops = dict(linewidth=2,
                color=[0.2, 0.2, 0.2])

ax.boxplot(y_data,
           positions=range(0, len(SNRCats)),
           showfliers=False,
           showcaps=False,
           medianprops=medianprops,
           whiskerprops=boxprops,
           boxprops=boxprops,
           widths=0.25)

# Add jittered dots ----------------------------------------------
xjitter = 0.03
yjitter = 0.06
x_data = [np.array([i] * len(d)) for i, d in enumerate(y_data)]
x_jittered = [x + stats.t(df=6, scale=xjitter).rvs(len(x)) for x in x_data]
y_jittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y)) for y in y_data]
for x, y, color in zip(x_jittered, y_jittered, mycolours[0:len(SNRCats)]):
    ax.scatter(x, y, s=10, color=color, alpha=0.2)

ax.set(yticks=range(0, 11), xticks=range(0, len(SNRCats)),
       xticklabels=SNRCats,
       xlabel=r"UAS $L_\mathrm{Aeq}$ - ambient $L_\mathrm{Aeq}$, dB",
       ylabel="Annoyance rating")
plt.show()

# segregated by ambient environment
fig, ax = plt.subplots(figsize=(8.15, 4.65))

data = partADataBySubj.loc[:, ['AmbientEnv', 'SNRlevel', 'Annoyance']]
data = data.sort_values('SNRlevel')

sns.violinplot(data=data, y='Annoyance', split=True, hue='AmbientEnv',
               x='SNRlevel', cut=0, palette=mycolours[0:2], inner='quart',
               width=0.5, bw_method='scott')

plt.setp(ax.collections, alpha=0.3)
# Add jittered dots ----------------------------------------------
xjitter = 0.03
yjitter = 0.06
y_subdataPark = list()
y_subdataStreet = list()
for ii in range(0, len(SNRCats)):
    y_subdataPark.append(partADataBySubj.loc[np.logical_and(
                                             partADataBySubj['AmbientEnv']
                                             == 'Park',
                                             partADataBySubj['SNRlevel']
                                             == SNRCats[ii]),
                                             'Annoyance'].values)
    y_subdataStreet.append(partADataBySubj.loc[np.logical_and(
                                               partADataBySubj['AmbientEnv']
                                               == 'Street',
                                               partADataBySubj['SNRlevel']
                                               == SNRCats[ii]),
                                               'Annoyance'].values)
x_subdataPark = [np.array([i] * len(d)) for i, d in enumerate(y_subdataPark)]
x_subdataStreet = [np.array([i] * len(d))
                   for i, d in enumerate(y_subdataStreet)]
x_Parkjittered = [x - 0.05 + stats.t(df=6, scale=xjitter).rvs(len(x))
                  for x in x_subdataPark]
y_Parkjittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y))
                  for y in y_subdataPark]
x_Streetjittered = [x + 0.05 + stats.t(df=6, scale=xjitter).rvs(len(x))
                    for x in x_subdataStreet]
y_Streetjittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y))
                    for y in y_subdataStreet]
for x, y in zip(x_Parkjittered, y_Parkjittered):
    ax.scatter(x, y, s=10, color=mycolours[0], alpha=0.2)
for x, y in zip(x_Streetjittered, y_Streetjittered):
    ax.scatter(x, y, s=10, color=mycolours[1], alpha=0.2)

ax.set(yticks=range(0, 11), xticks=range(0, len(SNRCats)), xticklabels=SNRCats,
       xlabel="UAS $L_\mathrm{Aeq}$ - ambient $L_\mathrm{Aeq}$, dB",
       ylabel="Annoyance rating", ylim=[-0.5, 10.5])
ax.legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=2, fontsize=11)
plt.show()

# by ops
fig, ax = plt.subplots(figsize=(5.25, 4))

opCats = ["No UAS", "Flyby", "Landing", "Takeoff"]

y_data = [partADataBySubj[partADataBySubj['UASOperation']
                          == op]['Annoyance'].values
          for op in opCats]

xjitter = 0.03
yjitter = 0.06
x_data = [np.array([i] * len(d)) for i, d in enumerate(y_data)]
x_jittered = [x + stats.t(df=6, scale=xjitter).rvs(len(x)) for x in x_data]
y_jittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y)) for y in y_data]

violins = ax.violinplot(y_data,
                        positions=range(0, len(opCats)),
                        widths=0.45,
                        bw_method='scott',
                        showmeans=False,
                        showmedians=False,
                        showextrema=False)
for ii, pc in enumerate(violins["bodies"]):
    pc.set_facecolor(mycolours[ii])
    pc.set_edgecolor([0.25, 0.25, 0.25])
    pc.set_linewidth(1)
    pc.set_alpha(0.25)

medianprops = dict(linewidth=4,
                   color=[0.2, 0.2, 0.2],
                   solid_capstyle="butt")

boxprops = dict(linewidth=2,
                color=[0.2, 0.2, 0.2])

ax.boxplot(y_data,
           positions=range(0, len(opCats)),
           showfliers=False,
           showcaps=False,
           medianprops=medianprops,
           whiskerprops=boxprops,
           boxprops=boxprops,
           widths=0.25)

# Add jittered dots ----------------------------------------------
for x, y, color in zip(x_jittered, y_jittered, mycolours[0:len(opCats)]):
    ax.scatter(x, y, s=10, color=color, alpha=0.2)

ax.set(yticks=range(0, 11), xticks=range(0, len(opCats)),
       xticklabels=opCats,
       xlabel=r"UAS operation",
       ylabel="Annoyance rating", ylim=[-0.5, 10.5])
plt.show()


# fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(8, 9.5))

# data = partADataBySubj.loc[:, ['UASOperation', 'AmbientEnv',
#                                'SNRlevel', 'Annoyance']]
# data = data.sort_values('SNRlevel')
# data0 = data.loc[data['UASOperation'].isin(['No UAS', 'Takeoff'])]
# data1 = data.loc[data['UASOperation'].isin(['No UAS', 'Flyby'])]
# data2 = data.loc[data['UASOperation'].isin(['No UAS', 'Landing'])]

# sns.violinplot(ax=axs[0], data=data0, y='Annoyance', split=True, hue='AmbientEnv',
#                x='SNRlevel', cut=0, palette=mycolours[0:2], inner='quart',
#                width=0.5, bw_method='scott', legend=True)
# sns.violinplot(ax=axs[1], data=data1, y='Annoyance', split=True, hue='AmbientEnv',
#                x='SNRlevel', cut=0, palette=mycolours[0:2], inner='quart',
#                width=0.5, bw_method='scott', legend=False)
# sns.violinplot(ax=axs[2], data=data2, y='Annoyance', split=True, hue='AmbientEnv',
#                x='SNRlevel', cut=0, palette=mycolours[0:2], inner='quart',
#                width=0.5, bw_method='scott', legend=False)

# plt.setp(axs[0].collections, alpha=0.3)
# plt.setp(axs[1].collections, alpha=0.3)
# plt.setp(axs[2].collections, alpha=0.3)

# axs[0].set(yticks=range(0, 11), xticks=range(0, len(SNRCats)),
#            xticklabels=SNRCats, xlabel=None, ylabel="Annoyance rating",
#            ylim=[-0.5, 10.5])
# axs[0].set_title("Takeoff", fontsize=11)
# axs[1].set(yticks=range(0, 11), xticks=range(0, len(SNRCats)),
#            xticklabels=SNRCats, xlabel=None, ylabel="Annoyance rating",
#            ylim=[-0.5, 10.5])
# axs[1].set_title("Flyby", fontsize=11)
# axs[2].set(yticks=range(0, 11), xticks=range(0, len(SNRCats)),
#            xticklabels=SNRCats,
#            xlabel="UAS $L_\mathrm{Aeq}$ - ambient $L_\mathrm{Aeq}$, dB",
#            ylabel="Annoyance rating", ylim=[-0.5, 10.5])
# axs[2].set_title("Landing", fontsize=11)
# axs[0].legend(bbox_to_anchor=(0.5, 1.3), loc='upper center', ncol=2,
#               fontsize=11)
# plt.show()

# by UAS type
fig, ax = plt.subplots(figsize=(5.25, 4))

vehicleCats = ["No UAS", "H520", "M300", "T150"]

y_data = [partADataBySubj[partADataBySubj['UASType']
                          == vehicle]['Annoyance'].values
          for vehicle in vehicleCats]

xjitter = 0.03
yjitter = 0.06
x_data = [np.array([i] * len(d)) for i, d in enumerate(y_data)]
x_jittered = [x + stats.t(df=6, scale=xjitter).rvs(len(x)) for x in x_data]
y_jittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y)) for y in y_data]

violins = ax.violinplot(y_data,
                        positions=range(0, len(vehicleCats)),
                        widths=0.45,
                        bw_method='scott',
                        showmeans=False,
                        showmedians=False,
                        showextrema=False)
for ii, pc in enumerate(violins["bodies"]):
    pc.set_facecolor(mycolours[ii])
    pc.set_edgecolor([0.25, 0.25, 0.25])
    pc.set_linewidth(1)
    pc.set_alpha(0.25)

medianprops = dict(linewidth=4,
                   color=[0.2, 0.2, 0.2],
                   solid_capstyle="butt")

boxprops = dict(linewidth=2,
                color=[0.2, 0.2, 0.2])

ax.boxplot(y_data,
           positions=range(0, len(vehicleCats)),
           showfliers=False,
           showcaps=False,
           medianprops=medianprops,
           whiskerprops=boxprops,
           boxprops=boxprops,
           widths=0.25)

# Add jittered dots ----------------------------------------------
for x, y, color in zip(x_jittered, y_jittered, mycolours[0:len(vehicleCats)]):
    ax.scatter(x, y, s=10, color=color, alpha=0.2)

ax.set(yticks=range(0, 11), xticks=range(0, len(vehicleCats)),
       xticklabels=vehicleCats,
       xlabel=r"UAS type",
       ylabel="Annoyance rating", ylim=[-0.5, 10.5])
plt.show()

# scatter plots
fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(11.5, 9.5))
data = dataByStimTestA.loc[dataByStimTestA.index.str.find("PwrScale")
                           == -1, ['UASOperation', 'AmbientEnv',
                                   'SNRlevel', 'ValenceMedian',
                                   'ArousalMedian', 'AnnoyMedian']]
data.sort_index(axis=0, inplace=True)

jitter_amount = 0.05
data_jittered = data.copy()
data_jittered[['ValenceMedian',
               'ArousalMedian',
               'AnnoyMedian']] = (data[['ValenceMedian',
                                       'ArousalMedian',
                                       'AnnoyMedian']]
                                       + np.random.normal(0,
                                                          jitter_amount,
                                                          data[['ValenceMedian',
                                                                'ArousalMedian',
                                                                'AnnoyMedian']].shape))

responseLabels = ["Median valence", "Median arousal", "Median annoyance"]
responseLims = [[1, 5], [1, 5], [0, 10]]
for ii in range(0, 3):
    for jj in range(0, 3):
        axs[ii, jj].plot([-1, 10], [-1, 10], color=[0.75, 0.75, 0.75],
                         alpha=0.5)

for ii, (response, colour) in enumerate(zip(['ArousalMedian', 'ValenceMedian',
                                             'AnnoyMedian'],
                                            [mycolours[0], mycolours[1],
                                             mycolours[2]])):
    for jj, operation in enumerate([['Flyby', 'Landing'], ['Flyby', 'Takeoff'],
                                    ['Takeoff', 'Landing']]):
        axs[ii, jj].scatter(data_jittered.loc[data['UASOperation'].isin(['No UAS',
                                                                         operation[0]])].loc[:,
                                                                                             response],
                            data_jittered.loc[data['UASOperation'].isin(['No UAS',
                                                                         operation[1]])].loc[:,
                                                                                             response],
                            alpha=0.2, color=colour)
        axs[ii, jj].set(xticks=range(0, 11), yticks=range(0, 11),
                        xlim = responseLims[ii], ylim=responseLims[ii],
                        xlabel=responseLabels[ii] + ", " + operation[0].lower(),
                        ylabel=responseLabels[ii] + ", " + operation[1].lower())

plt.show()

# detection bars



# Part B
# ------

# UAS LAeq

UASLAeqCats = list(partBDataBySubj['UASLAeq'].sort_values().unique())

# Annoyance

fig, ax = plt.subplots(figsize=(5, 4))

y_data = [partBDataBySubj[partBDataBySubj['UASLAeq']
                          == LAeq]['Annoyance'].values
          for LAeq in UASLAeqCats]

xjitter = 0.03
yjitter = 0.06
x_data = [np.array([i] * len(d)) for i, d in enumerate(y_data)]
x_jittered = [x + stats.t(df=6, scale=xjitter).rvs(len(x)) for x in x_data]
y_jittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y)) for y in y_data]

violins = ax.violinplot(y_data,
                        positions=range(0, len(UASLAeqCats)),
                        widths=0.45,
                        bw_method='scott',
                        showmeans=False,
                        showmedians=False,
                        showextrema=False)
for ii, pc in enumerate(violins["bodies"]):
    pc.set_facecolor(mycolours[ii])
    pc.set_edgecolor([0.25, 0.25, 0.25])
    pc.set_linewidth(1)
    pc.set_alpha(0.25)

medianprops = dict(linewidth=4,
                   color=[0.2, 0.2, 0.2],
                   solid_capstyle="butt")

boxprops = dict(linewidth=2,
                color=[0.2, 0.2, 0.2])

ax.boxplot(y_data,
           positions=range(0, len(UASLAeqCats)),
           showfliers=False,
           showcaps=False,
           medianprops=medianprops,
           whiskerprops=boxprops,
           boxprops=boxprops,
           widths=0.25)

# Add jittered dots ----------------------------------------------
for x, y, color in zip(x_jittered, y_jittered, mycolours[0:len(UASLAeqCats)]):
    ax.scatter(x, y, s=10, color=color, alpha=0.2)

ax.set(yticks=range(0, 11), xticks=range(0, len(UASLAeqCats)),
       xticklabels=UASLAeqCats,
       xlabel=r"UAS $L_\mathrm{Aeq}$, dB",
       ylabel="Annoyance rating", ylim=[-0.5, 10.5])
plt.show()

# UAS events

UASEventCats = list(partBDataBySubj['UASEvents'].sort_values().unique().astype(int))

# Annoyance

fig, ax = plt.subplots(figsize=(7, 4))

y_data = [partBDataBySubj[partBDataBySubj['UASEvents']
                          == UASEvents]['Annoyance'].values
          for UASEvents in UASEventCats]

xjitter = 0.03
yjitter = 0.06
x_data = [np.array([i] * len(d)) for i, d in enumerate(y_data)]
x_jittered = [x + stats.t(df=6, scale=xjitter).rvs(len(x)) for x in x_data]
y_jittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y)) for y in y_data]

violins = ax.violinplot(y_data,
                        positions=range(0, len(UASEventCats)),
                        widths=0.45,
                        bw_method='scott',
                        showmeans=False,
                        showmedians=False,
                        showextrema=False)
for ii, pc in enumerate(violins["bodies"]):
    pc.set_facecolor(mycolours[ii])
    pc.set_edgecolor([0.25, 0.25, 0.25])
    pc.set_linewidth(1)
    pc.set_alpha(0.25)

medianprops = dict(linewidth=4,
                   color=[0.2, 0.2, 0.2],
                   solid_capstyle="butt")

boxprops = dict(linewidth=2,
                color=[0.2, 0.2, 0.2])

ax.boxplot(y_data,
           positions=range(0, len(UASEventCats)),
           showfliers=False,
           showcaps=False,
           medianprops=medianprops,
           whiskerprops=boxprops,
           boxprops=boxprops,
           widths=0.25)

# Add jittered dots ----------------------------------------------
for x, y, color in zip(x_jittered, y_jittered, mycolours[0:len(UASEventCats)]):
    ax.scatter(x, y, s=10, color=color, alpha=0.2)

ax.set(yticks=range(0, 11), xticks=range(0, len(UASEventCats)),
       xticklabels=UASEventCats,
       xlabel=r"UAS events",
       ylabel="Annoyance rating", ylim=[-0.5, 10.5])
plt.show()

# segregated by LAeq

fig, ax = plt.subplots(figsize=(5.5, 4.65))

data = partBDataBySubj.loc[partBDataBySubj['UASEvents'] > 0,
                           ['UASEvents', 'UASLAeq', 'Annoyance']]
data = data.sort_values('UASEvents')
data['UASLAeq'] = pd.Categorical(data['UASLAeq'], ["54", "60"])

UASEventCats = list(data['UASEvents'].sort_values().unique().astype(int))
UASLAeqCats = list(data['UASLAeq'].sort_values().unique())

sns.violinplot(data=data, y='Annoyance', split=True, hue='UASLAeq',
               x='UASEvents', cut=0, palette=mycolours[0:len(UASLAeqCats)], inner='quart',
               width=0.5, bw_method='scott')

plt.setp(ax.collections, alpha=0.3)
# Add jittered dots ----------------------------------------------
y_subdata54 = list()
y_subdata60 = list()
for ii in range(0, len(UASEventCats)):
    y_subdata54.append(partBDataBySubj.loc[np.logical_and(
                                           partBDataBySubj['UASLAeq']
                                           == "54",
                                           partBDataBySubj['UASEvents']
                                           == UASEventCats[ii]),
                       'Annoyance'].values)
    y_subdata60.append(partBDataBySubj.loc[np.logical_and(
                                           partBDataBySubj['UASLAeq']
                                           == "60",
                                           partBDataBySubj['UASEvents']
                                           == UASEventCats[ii]),
                       'Annoyance'].values)
x_subdata54 = [np.array([i] * len(d)) for i, d in enumerate(y_subdata54)]
x_subdata60 = [np.array([i] * len(d))
                   for i, d in enumerate(y_subdata60)]
x_54jittered = [x - 0.05 + stats.t(df=6, scale=xjitter).rvs(len(x))
                  for x in x_subdata54]
y_54jittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y))
                  for y in y_subdata54]
x_60jittered = [x + 0.05 + stats.t(df=6, scale=xjitter).rvs(len(x))
                    for x in x_subdata60]
y_60jittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y))
                    for y in y_subdata60]
for x, y in zip(x_54jittered, y_54jittered):
    ax.scatter(x, y, s=10, color=mycolours[0], alpha=0.2)
for x, y in zip(x_60jittered, y_60jittered):
    ax.scatter(x, y, s=10, color=mycolours[1], alpha=0.2)

ax.set(yticks=range(0, 11), xticks=range(0, len(UASEventCats)),
       xticklabels=UASEventCats, xlabel=r"UAS Events",
       ylabel="Annoyance rating", ylim=[-0.5, 10.5])
ax.legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=2, fontsize=11,
          title=r"UAS $L_\mathrm{Aeq}$, dB")
plt.show()

# plot participant sample feature data

# age distribution
fig, ax = plt.subplots(figsize=(4, 3))

plt.hist(postTestResponses.Age[postTestResponses.Age != "No answer"],
         bins=np.arange(18.5, 62.5, 4), histtype='step', color=mycolours[0])
ax.set(xticks=range(18, 62, 2), xlim=[18, 60], yticks=range(0, 12),
       xlabel="Age, years", ylabel="No. of participants")
plt.show()

# residential area
fig, ax = plt.subplots(figsize=(4, 3))

areaCats = np.sort(postTestResponses['Home Area'].unique())
data = postTestResponses[['Home Area']
                         + ['Area soundscape']].apply('.'.join,
                                                      axis=1).value_counts().sort_index()
scapeCounts = {"Calm": np.array([1, 8, 11]),
               "Chaotic": np.array([0, 0, 4]),
               "Monotonous": np.array([0, 2, 2]),
               "Vibrant": np.array([0, 4, 10])}
# wedges, texts = ax.pie(postTestResponses['Home Area'].value_counts(),
#                        wedgeprops=dict(width=size, edgecolor='w'),
#                        startangle=-40, radius=1, colors=inner_colors)

width = 0.5

bottom = np.zeros(3)

for ii, (boolean, scapeCounts) in enumerate(scapeCounts.items()):
    p = ax.bar(areaCats, scapeCounts,
               label=boolean, bottom=bottom, color=mycolours[ii])
    bottom += scapeCounts

ax.legend(loc="upper left", title="Area soundscape", fontsize=9)
ax.set(yticks=range(0, 32, 2), ylabel="No. of participants", xlabel="Area of residence")

plt.show()
