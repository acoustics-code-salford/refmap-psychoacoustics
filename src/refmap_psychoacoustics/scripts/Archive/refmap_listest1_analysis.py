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
import pingouin as pg

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
             (251, 164, 49), (204, 153, 255), (90, 192, 255), (80, 245, 233),
             (255, 90, 192), (164, 201, 242), (255, 254, 139), (255, 243, 255)]
mycolours = [tuple(shade/255 for shade in colour) for colour in mycolours]


# enable copy-on-write mode for Pandas (will be default from Pandas 3.0)
pd.options.mode.copy_on_write = True

# import data
app = QApplication(sys.argv)
fileExts = "*.csv"

# data by stimulus
dataByStimFilePath = list(QFileDialog.getOpenFileName(filter=fileExts,
                                                      caption=r"Open refmap_listest1_testdata_ByStim.csv in: \03 Experiment\Experiment 1\Analysis\PostProcess"))[0]
dataByStimTest = pd.read_csv(dataByStimFilePath, index_col=0)

# Part A
dataByStimAFilePath = list(QFileDialog.getOpenFileName(filter=fileExts,
                                                       caption=r"Open refmap_listest1_testdataA_ByStim.csv in: \03 Experiment\Experiment 1\Analysis\PostProcess"))[0]
dataByStimTestA = pd.read_csv(dataByStimAFilePath, index_col=0)

# Part A notice data subselection
dataByStimANoticeFilePath = list(QFileDialog.getOpenFileName(filter=fileExts,
                                                             caption=r"Open refmap_listest1_testdataANoticeFilt_ByStim.csv in: \03 Experiment\Experiment 1\Analysis\PostProcess"))[0]
dataByStimTestANotice = pd.read_csv(dataByStimANoticeFilePath, index_col=0)

# Part B
dataByStimBFilePath = list(QFileDialog.getOpenFileName(filter=fileExts,
                                                       caption=r"Open refmap_listest1_testdataB_ByStim.csv in: \03 Experiment\Experiment 1\Analysis\PostProcess"))[0]
dataByStimTestB = pd.read_csv(dataByStimBFilePath, index_col=0)

# data by subject
allDataBySubjFilePath = list(QFileDialog.getOpenFileName(filter=fileExts,
                                                         caption=r"Open refmap_listest1_testdata_BySubj.csv in: \03 Experiment\Experiment 1\Analysis\PostProcess"))[0]
allDataBySubj = pd.read_csv(allDataBySubjFilePath, index_col=0)

# Part A
partADataFilePath = list(QFileDialog.getOpenFileName(filter=fileExts,
                                                     caption=r"Open refmap_listest1_testdataA_BySubj.csv in: \03 Experiment\Experiment 1\Analysis\PostProcess"))[0]
partADataBySubj = pd.read_csv(partADataFilePath, index_col=False)

# Part A notice data subselection
partANoticeDataFilePath = list(QFileDialog.getOpenFileName(filter=fileExts,
                                                           caption=r"Open refmap_listest1_testdataANoticeFilt_BySubj.csv in: \03 Experiment\Experiment 1\Analysis\PostProcess"))[0]
partANoticeDataBySubj = pd.read_csv(partANoticeDataFilePath, index_col=False)

# Part B
partBDataFilePath = list(QFileDialog.getOpenFileName(filter=fileExts,
                                                     caption=r"Open refmap_listest1_testdataB_BySubj.csv in: \03 Experiment\Experiment 1\Analysis\PostProcess"))[0]
partBDataBySubj = pd.read_csv(partBDataFilePath, index_col=False)

# pre and post test data
preTestDataFilePath = list(QFileDialog.getOpenFileName(filter=fileExts,
                                                       caption=r"Open refmap_listest1_pretestdata.csv in: 03 Experiment\Experiment 1\Analysis\PostProcess"))[0]
preTestResponses = pd.read_csv(preTestDataFilePath, index_col=0)
postTestDataFilePath = list(QFileDialog.getOpenFileName(filter=fileExts,
                                                        caption=r"Open refmap_listest1_posttestdata.csv in: 03 Experiment\Experiment 1\Analysis\PostProcess"))[0]
postTestResponses = pd.read_csv(postTestDataFilePath, index_col=0)

# categorise columns
SNRCats = ["No UAS", "-16", "-10", "-4", "2", "8"]
UASLAeqCats = ["No UAS", "42", "48", "54", "60"]
opCats = ["No UAS", "Flyby", "Landing", "Takeoff"]
vehicleCats = ["No UAS", "H520", "M300", "T150"]

for dataset in [dataByStimTest, dataByStimTestA, dataByStimTestANotice,
                allDataBySubj, partADataBySubj, partANoticeDataBySubj]:
    dataset['SNRlevel'] = pd.Categorical(dataset['SNRlevel'], SNRCats)
    dataset['UASLAeq'] = pd.Categorical(dataset['UASLAeq'], UASLAeqCats)
    dataset['UASOperation'] = pd.Categorical(dataset['UASOperation'], opCats)
    dataset['UASType'] = pd.Categorical(dataset['UASType'], vehicleCats)

SNRCatsB = ["No UAS", "2", "8"]
UASLAeqCatsB = ["No UAS", "54", "60"]
opCatsB = ["No UAS", "Flyby"]
vehicleCatsB = ["No UAS", "H520", "T150"]

for dataset in [dataByStimTestB, partBDataBySubj]:
    dataset['SNRlevel'] = pd.Categorical(dataset['SNRlevel'], SNRCatsB)
    dataset['UASLAeq'] = pd.Categorical(dataset['UASLAeq'], UASLAeqCatsB)
    dataset['UASOperation'] = pd.Categorical(dataset['UASOperation'], opCatsB)
    dataset['UASType'] = pd.Categorical(dataset['UASType'], vehicleCatsB)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# EXPLORATORY ANALYSIS
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# -----------------------
# Anomalous data checking
# -----------------------

# Tukey outlier test
kFact = 1.5  # multiplying factor to identify outliers from IQ range

# Annoyance
AnnoyLowQ = dataByStimTest.loc[:,
                               'Annoyance_1':'Annoyance_19'].quantile(0.25,
                                                                      axis=1).to_frame()
AnnoyHiQ = dataByStimTest.loc[:,
                              'Annoyance_1':'Annoyance_19'].quantile(0.75,
                                                                     axis=1).to_frame()
AnnoyIQ = AnnoyHiQ.values - AnnoyLowQ.values

AnnoyLowOut = AnnoyLowQ.values - kFact*(AnnoyIQ)
AnnoyHiOut = AnnoyHiQ.values + kFact*(AnnoyIQ)

AnnoyOutTest = pd.DataFrame(data=np.logical_or(dataByStimTest.loc[:,
                                                                  'Annoyance_1':'Annoyance_19'].values
                                               < AnnoyLowOut,
                                               dataByStimTest.loc[:,
                                                                  'Annoyance_1':'Annoyance_19'].values
                                               > AnnoyHiOut),
                            index=dataByStimTest.loc[:,
                                                     'Annoyance_1':'Annoyance_19'].index,
                            columns=[participant.replace("Annoyance_", "")
                                     for participant
                                     in dataByStimTest.loc[:,
                                                           'Annoyance_1':'Annoyance_19'].columns])
AnnoyOutTestScore = AnnoyOutTest.sum(axis=0)
AnnoyOutTestScore.index = AnnoyOutTestScore.index.astype(dtype=int,
                                                         copy=False)
AnnoyOutTestScore.sort_index(inplace=True)

# Valence
ValenceLowQ = dataByStimTest.loc[:,
                                 'Valence_1':'Valence_19'].quantile(0.25,
                                                                    axis=1).to_frame()
ValenceHiQ = dataByStimTest.loc[:,
                                'Valence_1':'Valence_19'].quantile(0.75,
                                                                   axis=1).to_frame()
ValenceIQ = ValenceHiQ.values - ValenceLowQ.values

ValenceLowOut = ValenceLowQ.values - kFact*(ValenceIQ)
ValenceHiOut = ValenceHiQ.values + kFact*(ValenceIQ)

ValenceOutTest = pd.DataFrame(data=np.logical_or(dataByStimTest.loc[:,
                                                                    'Valence_1':'Valence_19'].values
                                                 < ValenceLowOut,
                                               dataByStimTest.loc[:,
                                                                  'Valence_1':'Valence_19'].values
                                               > ValenceHiOut),
                            index=dataByStimTest.loc[:,
                                                     'Valence_1':'Valence_19'].index,
                            columns=[participant.replace("Valence_", "") for participant
                                     in dataByStimTest.loc[:,
                                                           'Valence_1':'Valence_19'].columns])
ValenceOutTestScore = ValenceOutTest.sum(axis=0)
ValenceOutTestScore.index = ValenceOutTestScore.index.astype(dtype=int,
                                                             copy=False)
ValenceOutTestScore.sort_index(inplace=True)

# Arousal
ArousalLowQ = dataByStimTest.loc[:,
                               'Arousal_1':'Arousal_19'].quantile(0.25,
                                                                      axis=1).to_frame()
ArousalHiQ = dataByStimTest.loc[:,
                              'Arousal_1':'Arousal_19'].quantile(0.75,
                                                                     axis=1).to_frame()
ArousalIQ = ArousalHiQ.values - ArousalLowQ.values

ArousalLowOut = ArousalLowQ.values - kFact*(ArousalIQ)
ArousalHiOut = ArousalHiQ.values + kFact*(ArousalIQ)

ArousalOutTest = pd.DataFrame(data=np.logical_or(dataByStimTest.loc[:,
                                                                    'Arousal_1':'Arousal_19'].values
                                                 < ArousalLowOut,
                                               dataByStimTest.loc[:,
                                                                  'Arousal_1':'Arousal_19'].values
                                               > ArousalHiOut),
                              index=dataByStimTest.loc[:,
                                                     'Arousal_1':'Arousal_19'].index,
                              columns=[participant.replace("Arousal_", "") for participant
                                       in dataByStimTest.loc[:,
                                                             'Arousal_1':'Arousal_19'].columns])
ArousalOutTestScore = ArousalOutTest.sum(axis=0)
ArousalOutTestScore.index = ArousalOutTestScore.index.astype(dtype=int,
                                                             copy=False)
ArousalOutTestScore.sort_index(inplace=True)

# Total outliers per participant
TotalOutTestScore = AnnoyOutTestScore + ValenceOutTestScore + ArousalOutTestScore

# plot response outliers
fig, ax = plt.subplots(figsize=(13, 4))
OutScoreProps = {"Annoyance": 100*AnnoyOutTestScore/AnnoyOutTest.shape[0],
                 "Valence": 100*ValenceOutTestScore/ValenceOutTest.shape[0],
                 "Arousal": 100*ArousalOutTestScore/ArousalOutTest.shape[0]}
width = 0.5
bottom = np.zeros(ArousalOutTestScore.index.size)

for ii, (responseType, OutScoreProp) in enumerate(OutScoreProps.items()):
    p = ax.bar(OutScoreProp.index.astype(str), OutScoreProp, width,
               label=responseType,
               bottom=bottom, color=mycolours[ii])
    bottom += OutScoreProp

barLabels = TotalOutTestScore.astype(str)
barLabels[barLabels != "0"] = ""
ax.bar_label(container=p, labels=barLabels, padding =3)
ax.legend(title="Response")
ax.set(xticks=OutScoreProp.index.astype(str), yticks=range(0, 110, 10), ylabel="Proportion of response outliers, %",
       xlabel="Participant ID#")
plt.show()

# Noticeability
subdata = partADataBySubj.loc[:, ['ID#', 'UAS_noticed',
                                  'UASPartLoudGMSTPowAvg']]
allIDs = subdata['ID#'].unique()
noticeCheck = pd.DataFrame(index=allIDs, columns=["loCount", "loTot",
                                                  "hiCount", "hiTot"])
for ID in allIDs:
    noticeCheck.loc[ID, 'loCount'] = subdata.loc[subdata['ID#']
                                                 == ID].loc[subdata['UASPartLoudGMSTPowAvg']
                                                            < 0.3].loc[:,
                                                                        'UAS_noticed'].sum()
    noticeCheck.loc[ID, 'loTot'] = np.size(subdata.loc[subdata['ID#']
                                                       == ID].loc[subdata['UASPartLoudGMSTPowAvg']
                                                                  < 0.3].loc[:, 'UAS_noticed'])
    noticeCheck.loc[ID, 'hiCount'] = (subdata.loc[subdata['ID#']
                                                  == ID].loc[subdata['UASPartLoudGMSTPowAvg']
                                                             > 13].loc[:, 'UAS_noticed'] == 0).sum()
    noticeCheck.loc[ID, 'hiTot'] = np.size(subdata.loc[subdata['ID#']
                                                       == ID].loc[subdata['UASPartLoudGMSTPowAvg']
                                                                  > 13].loc[:, 'UAS_noticed'])

noticeCheck['lohiCount'] = noticeCheck['loCount'] + noticeCheck['hiCount']
noticeCheck['lohiTot'] = noticeCheck['loTot'] + noticeCheck['hiTot']
noticeCheck['loProp'] = 100*noticeCheck['loCount'] / noticeCheck['loTot']
noticeCheck['hiProp'] = 100*noticeCheck['hiCount'] / noticeCheck['hiTot']

fig, ax = plt.subplots(figsize=(13, 4))

x = np.arange(len(noticeCheck.index))  # the label locations
width = 0.45  # the width of the bars
multiplier = 0

labels = ["$\\Delta N < 0.3$ sones", "$\\Delta N > 13$ sones"]
for ii, lohi in enumerate(noticeCheck[['loProp', 'hiProp']].columns):
    offset = width * multiplier
    rects = ax.bar(x + offset, noticeCheck[lohi], width, label=labels[ii],
                   color=mycolours[ii])
    barLabels = np.array([round(value) for value in rects.datavalues], dtype=str)
    barLabels[barLabels != "0"] = ""
    ax.bar_label(container=rects, labels=barLabels, padding=3)
    multiplier += 1

ax.set(xticks=np.arange(len(noticeCheck.index)) + width/2,
       xticklabels=noticeCheck.index.astype(str), yticks=range(0, 110, 10),
       ylabel="Proportion of suspect 'UAS noticed' responses, %",
       xlabel="Participant ID#")
ax.legend()
plt.show()


# ------------------------------------
# plot participant sample feature data
# ------------------------------------

# age distribution
fig, ax = plt.subplots(figsize=(4, 3))

plt.hist(x=postTestResponses.Age[postTestResponses.Age != "No answer"].astype(int),
         bins=np.arange(18.5, 62.5, 4), histtype='step', color=mycolours[0])
ax.set(xticks=range(18, 62, 2), xlim=[18, 60], yticks=range(0, 12),
       xlabel="Age, years", ylabel="No. of participants")
plt.show()

# residential area
fig, ax = plt.subplots(figsize=(4, 3))

areaCats = np.sort(postTestResponses['Home_Area'].unique())
data = postTestResponses[['Home_Area']
                         + ['Area_soundscape']].apply('.'.join,
                                                      axis=1).value_counts().sort_index()
scapeCounts = {"Calm": np.array([1, 8, 11]),
               "Chaotic": np.array([0, 0, 4]),
               "Monotonous": np.array([0, 2, 2]),
               "Vibrant": np.array([0, 4, 10])}

width = 0.5

bottom = np.zeros(3)

for ii, (boolean, scapeCounts) in enumerate(scapeCounts.items()):
    p = ax.bar(areaCats, scapeCounts,
               label=boolean, bottom=bottom, color=mycolours[ii])
    bottom += scapeCounts

ax.legend(loc="upper left", title="Area soundscape", fontsize=9)
ax.set(yticks=range(0, 32, 2), ylabel="No. of participants",
       xlabel="Area of residence")

plt.show()

# noise sensitivity distribution
fig, ax = plt.subplots(figsize=(5, 3))

NSSnorm = 5  # normalisation constant to reduce all scores to 0-5 scale

plt.hist(x=postTestResponses.NSSTotal - NSSnorm, bins=np.arange(4.5, 30.5, 1),
         histtype='step', color=mycolours[1])
ax.set(xticks=range(0, 26, 1), xlim=[-0.5, 25.5], yticks=range(0, 6),
       xlabel="Noise sensitivity score (0-25)", ylabel="No. of participants")
plt.show()

# AAM attitude

labels = postTestResponses['AAM_attitude'].value_counts().index
sizes = postTestResponses['AAM_attitude'].value_counts().values

fig, ax = plt.subplots(figsize=(4, 4))
ax.pie(sizes, labels=labels, wedgeprops=dict(width=0.5), startangle=-90,
       colors=mycolours)
plt.show()

# -----------------------------------------------------------------------------
# Exploratory plots: All data
# -----------------------------------------------------------------------------

# median vs mean annoyance
fig, ax = plt.subplots(figsize=(4.5, 4))

ax.plot([0, 10], [0, 10], linestyle='-', color=[0.25, 0.25, 0.25], alpha=0.25)
ax.scatter(x=dataByStimTestA['AnnoyMedian'], y=dataByStimTestA['AnnoyMean'],
           color=mycolours[2], alpha=0.25, label="Part A", marker='x')
ax.set(xticks=range(0, 11), xlabel="Median annoyance rating",
       yticks=range(0, 11), ylabel="Mean annoyance rating")

# linear regression line
b, a = np.polyfit(dataByStimTestA['AnnoyMedian'], dataByStimTestA['AnnoyMean'],
                  deg=1)
xseq = np.linspace(0, 10, num=11)
ax.plot(xseq, a + b * xseq, color=mycolours[2], linewidth=1)

ax.scatter(x=dataByStimTestB['AnnoyMedian'], y=dataByStimTestB['AnnoyMean'],
           color=mycolours[5], alpha=0.35, label="Part B", marker='o',
           facecolor='none')
ax.set(xticks=range(0, 11), xlabel="Median annoyance rating",
       yticks=range(0, 11), ylabel="Mean annoyance rating")

# linear regression line
b, a = np.polyfit(dataByStimTestB['AnnoyMedian'], dataByStimTestB['AnnoyMean'],
                  deg=1)
xseq = np.linspace(0, 10, num=11)
ax.plot(xseq, a + b * xseq, color=mycolours[5], linewidth=1)


ax.grid(alpha=0.15, linestyle='--')
ax.legend()
plt.show()

# histogram plot of extreme annoyance distributions
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(11, 4))

data = [partADataBySubj.loc[partADataBySubj['Recording']
                            == "A2_CALBIN_Pa.wav", 'Annoyance'],
        partADataBySubj.loc[partADataBySubj['Recording']
                            == "A2_H520_F_2_R_CALBIN_Pa.wav", 'Annoyance'],
        partADataBySubj.loc[partADataBySubj['Recording']
                            == "A2_H520_L_1_CALBIN_Pa.wav", 'Annoyance']]

for ii, ax in enumerate(axs):
    sns.histplot(data=data[ii], color=mycolours[0], alpha=0.25,
                 bins=np.arange(-0.5, 11.5, 1), ax=ax, kde=True)
    ax.vlines(x=np.median(data[ii]), ymin=0, ymax=20,
              colors=mycolours[1], linewidth=1.5)
    ax.vlines(x=np.mean(data[ii]), ymin=0, ymax=20,
              colors=mycolours[2], linewidth=1.5)
    ax.set(xticks=range(0, 11, 1), ylim=[0, 16], xlabel="Annoyance rating")
    ax.legend(labels=["Data KDE", "Median", "Mean", "Data"])
    if ii == 0:
        ax.set(ylabel="Response count")
    else:
        ax.set(ylabel=None)

# TODO
# loudness index comparisons


# -----------------------------------------------------------------------------
# Exploratory plots: Part A
# -----------------------------------------------------------------------------

# Response analysis: Annoyance
# ----------------------------

# UAS LAeq

UASLAeqCats = list(partADataBySubj['UASLAeq'].sort_values().unique())

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

# Add jittered dots
for x, y, color in zip(x_jittered, y_jittered, mycolours[0:len(UASLAeqCats)]):
    ax.scatter(x, y, s=10, color=color, alpha=0.2)

ax.set(yticks=range(0, 11), xticks=range(0, len(UASLAeqCats)),
       xticklabels=UASLAeqCats,
       xlabel="UAS $L_\\mathrm{Aeq}$, dB",
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
# Add jittered dots
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
       xticklabels=UASLAeqCats, xlabel="UAS $L_\\mathrm{Aeq}$, dB",
       ylabel="Annoyance rating", ylim=[-0.5, 10.5])
ax.legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=2, fontsize=11)
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

# Add jittered dots
xjitter = 0.03
yjitter = 0.06
x_data = [np.array([i] * len(d)) for i, d in enumerate(y_data)]
x_jittered = [x + stats.t(df=6, scale=xjitter).rvs(len(x)) for x in x_data]
y_jittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y)) for y in y_data]
for x, y, color in zip(x_jittered, y_jittered, mycolours[0:len(SNRCats)]):
    ax.scatter(x, y, s=10, color=color, alpha=0.2)

ax.set(yticks=range(0, 11), xticks=range(0, len(SNRCats)),
       xticklabels=SNRCats,
       xlabel="UAS $L_\\mathrm{Aeq}$ - ambient $L_\\mathrm{Aeq}$, dB",
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
# Add jittered dots
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
       xlabel="UAS $L_\\mathrm{Aeq}$ - ambient $L_\\mathrm{Aeq}$, dB",
       ylabel="Annoyance rating", ylim=[-0.5, 10.5])
ax.legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=2, fontsize=11)
plt.show()

# Flight ops
# ----------
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

# Add jittered dots
for x, y, color in zip(x_jittered, y_jittered, mycolours[0:len(opCats)]):
    ax.scatter(x, y, s=10, color=color, alpha=0.2)

ax.set(yticks=range(0, 11), xticks=range(0, len(opCats)),
       xticklabels=opCats,
       xlabel="UAS operation",
       ylabel="Annoyance rating", ylim=[-0.5, 10.5])
plt.show()


# fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(8, 9.5))

# data = partADataBySubj.loc[:, ['UASOperation', 'AmbientEnv',
#                                'SNRlevel', 'Annoyance']]
# data = data.sort_values('SNRlevel')
# data0 = data.loc[data['UASOperation'].isin(['No UAS', 'Takeoff'])]
# data1 = data.loc[data['UASOperation'].isin(['No UAS', 'Flyby'])]
# data2 = data.loc[data['UASOperation'].isin(['No UAS', 'Landing'])]

# sns.violinplot(ax=axs[0], data=data0, y='Annoyance', split=True,
#                hue='AmbientEnv', x='SNRlevel', cut=0,
#                palette=mycolours[0:2], inner='quart',
#                width=0.5, bw_method='scott', legend=True)
# sns.violinplot(ax=axs[1], data=data1, y='Annoyance', split=True,
#                hue='AmbientEnv', x='SNRlevel', cut=0,
#                palette=mycolours[0:2], inner='quart',
#                width=0.5, bw_method='scott', legend=False)
# sns.violinplot(ax=axs[2], data=data2, y='Annoyance', split=True,
#                hue='AmbientEnv', x='SNRlevel', cut=0,
#                palette=mycolours[0:2], inner='quart',
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
#            xlabel="UAS $L_\\mathrm{Aeq}$ - ambient $L_\\mathrm{Aeq}$, dB",
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

# Add jittered dots
for x, y, color in zip(x_jittered, y_jittered, mycolours[0:len(vehicleCats)]):
    ax.scatter(x, y, s=10, color=color, alpha=0.2)

ax.set(yticks=range(0, 11), xticks=range(0, len(vehicleCats)),
       xticklabels=vehicleCats,
       xlabel="UAS type",
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
                        xlim=responseLims[ii], ylim=responseLims[ii],
                        xlabel=responseLabels[ii] + ", "
                        + operation[0].lower(),
                        ylabel=responseLabels[ii] + ", "
                        + operation[1].lower())

plt.show()


# Response analysis: Valence & arousal
# ------------------------------------

# UAS LAeq

UASLAeqCats = list(partADataBySubj['UASLAeq'].sort_values().unique())

# segregated by ambient environment
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7, 7))

data0 = partADataBySubj.loc[:, ['AmbientEnv', 'UASLAeq', 'Valence']]
data0 = data0.sort_values('UASLAeq')
data1 = partADataBySubj.loc[:, ['AmbientEnv', 'UASLAeq', 'Arousal']]
data1 = data1.sort_values('UASLAeq')

sns.violinplot(ax=axs[0], data=data0, y='Valence', split=True,
               hue='AmbientEnv', x='UASLAeq', cut=0,
               palette=mycolours[0:2], inner='quart',
               width=0.5, bw_method=0.55, legend=True)
sns.violinplot(ax=axs[1], data=data1, y='Arousal', split=True,
               hue='AmbientEnv', x='UASLAeq', cut=0,
               palette=mycolours[0:2], inner='quart',
               width=0.5, bw_method=0.55, legend=False)

plt.setp(axs[0].collections, alpha=0.3)
plt.setp(axs[1].collections, alpha=0.3)
# Add jittered dots
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
           xticklabels=UASLAeqCats, xlabel="UAS $L_\\mathrm{Aeq}$, dB",
           ylabel="Arousal rating", ylim=[0.5, 5.5])
axs[0].legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=2,
              fontsize=11)
plt.show()


# Perception analysis: Noticeability
# ----------------------------------

# unfiltered participant selection data
fig, ax = plt.subplots(figsize=(7.15, 4))

SNRCats = list(dataByStimTestA['SNRlevel'].sort_values().unique())

y_data = [dataByStimTestA.loc[dataByStimTestA['SNRlevel']
                              == SNR, 'NoticedProportion'].values
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

# Add jittered dots
xjitter = 0.02
yjitter = 0.0
x_data = [np.array([i] * len(d)) for i, d in enumerate(y_data)]
x_jittered = [x + stats.t(df=6, scale=xjitter).rvs(len(x)) for x in x_data]
y_jittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y)) for y in y_data]
for x, y, color in zip(x_jittered, y_jittered, mycolours[0:len(SNRCats)]):
    ax.scatter(x, y, s=10, color=color, alpha=0.2)

ax.set(yticks=np.arange(0, 1.1, 0.1), yticklabels=range(0, 110, 10),
       xticks=range(0, len(SNRCats)), xticklabels=SNRCats,
       xlabel="UAS $L_\\mathrm{Aeq}$ - ambient $L_\\mathrm{Aeq}$, dB",
       ylabel="Proportion UAS noticed, %")
plt.show()


# segregated by ambient environment
fig, ax = plt.subplots(figsize=(7.15, 4.65))

data = dataByStimTestA.loc[:, ['AmbientEnv', 'SNRlevel', 'NoticedProportion']]
data = data.sort_values('AmbientEnv', ignore_index=True)
data = data.sort_values('SNRlevel')


sns.violinplot(data=data, y='NoticedProportion', split=True, hue='AmbientEnv',
               x='SNRlevel', cut=0, palette=mycolours[0:2], inner='quart',
               width=0.65, bw_method='scott')

plt.setp(ax.collections, alpha=0.3)
# Add jittered dots
xjitter = 0.01
yjitter = 0.0
y_subdataPark = list()
y_subdataStreet = list()
for ii in range(0, len(SNRCats)):
    y_subdataPark.append(dataByStimTestA.loc[np.logical_and(
                                             dataByStimTestA['AmbientEnv']
                                             == 'Park',
                                             dataByStimTestA['SNRlevel']
                                             == SNRCats[ii]),
                                             'NoticedProportion'].values)
    y_subdataStreet.append(dataByStimTestA.loc[np.logical_and(
                                               dataByStimTestA['AmbientEnv']
                                               == 'Street',
                                               dataByStimTestA['SNRlevel']
                                               == SNRCats[ii]),
                                               'NoticedProportion'].values)
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
    ax.scatter(x, y, s=10, color=mycolours[0], alpha=0.3)
for x, y in zip(x_Streetjittered, y_Streetjittered):
    ax.scatter(x, y, s=10, color=mycolours[1], alpha=0.3)

ax.set(yticks=np.arange(0, 1.1, 0.1), yticklabels=range(0, 110, 10),
       xticks=range(0, len(SNRCats)), xticklabels=SNRCats,
       xlabel="UAS $L_\\mathrm{Aeq}$ - ambient $L_\\mathrm{Aeq}$, dB",
       ylabel="Proportion UAS noticed, %", ylim=[-0.05, 1.05])
ax.legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=2, fontsize=11)
plt.show()

# re-run analysis with recalculated proportion from selected participants
fig, ax = plt.subplots(figsize=(7.15, 4))

SNRCats = list(dataByStimTestA['SNRlevel'].sort_values().unique())

omitParticipants = [2, 5, 21, 32, 34, 36, 38, 39, 40, 44, 45]
omitColumns = ["UAS_noticed_" + str(partID) for partID in omitParticipants]
keepColumns = [col for col in dataByStimTestA.columns[dataByStimTestA.columns.str.find("UAS_noticed_") == 0]
               if col not in omitColumns]
dataByStimTestA['NoticedTotalFilt'] = dataByStimTestA[keepColumns].sum(axis=1)
dataByStimTestA['NoticedPropFilt'] = dataByStimTestA[keepColumns].sum(axis=1)/len(keepColumns)

y_data = [dataByStimTestA.loc[dataByStimTestA['SNRlevel']
                              == SNR, 'NoticedPropFilt'].values
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

# Add jittered dots
xjitter = 0.02
yjitter = 0.0
x_data = [np.array([i] * len(d)) for i, d in enumerate(y_data)]
x_jittered = [x + stats.t(df=6, scale=xjitter).rvs(len(x)) for x in x_data]
y_jittered = [y + stats.t(df=6, scale=yjitter).rvs(len(y)) for y in y_data]
for x, y, color in zip(x_jittered, y_jittered, mycolours[0:len(SNRCats)]):
    ax.scatter(x, y, s=10, color=color, alpha=0.2)

ax.set(yticks=np.arange(0, 1.1, 0.1), yticklabels=range(0, 110, 10),
       xticks=range(0, len(SNRCats)), xticklabels=SNRCats,
       xlabel="UAS $L_\\mathrm{Aeq}$ - ambient $L_\\mathrm{Aeq}$, dB",
       ylabel="Proportion UAS noticed, %")
plt.show()


# segregated by ambient environment
fig, ax = plt.subplots(figsize=(7.15, 4.65))

data = dataByStimTestA.loc[:, ['AmbientEnv', 'SNRlevel', 'NoticedPropFilt']]
data = data.sort_values('AmbientEnv', ignore_index=True)
data = data.sort_values('SNRlevel')


sns.violinplot(data=data, y='NoticedPropFilt', split=True, hue='AmbientEnv',
               x='SNRlevel', cut=0, palette=mycolours[0:2], inner='quart',
               width=0.65, bw_method='scott')

plt.setp(ax.collections, alpha=0.3)
# Add jittered dots
xjitter = 0.01
yjitter = 0.0
y_subdataPark = list()
y_subdataStreet = list()
for ii in range(0, len(SNRCats)):
    y_subdataPark.append(dataByStimTestA.loc[np.logical_and(
                                             dataByStimTestA['AmbientEnv']
                                             == 'Park',
                                             dataByStimTestA['SNRlevel']
                                             == SNRCats[ii]),
                                             'NoticedPropFilt'].values)
    y_subdataStreet.append(dataByStimTestA.loc[np.logical_and(
                                               dataByStimTestA['AmbientEnv']
                                               == 'Street',
                                               dataByStimTestA['SNRlevel']
                                               == SNRCats[ii]),
                                               'NoticedPropFilt'].values)
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
    ax.scatter(x, y, s=10, color=mycolours[0], alpha=0.3)
for x, y in zip(x_Streetjittered, y_Streetjittered):
    ax.scatter(x, y, s=10, color=mycolours[1], alpha=0.3)

ax.set(yticks=np.arange(0, 1.1, 0.1), yticklabels=range(0, 110, 10),
       xticks=range(0, len(SNRCats)), xticklabels=SNRCats,
       xlabel="UAS $L_\\mathrm{Aeq}$ - ambient $L_\\mathrm{Aeq}$, dB",
       ylabel="Proportion UAS noticed, %", ylim=[-0.05, 1.05])
ax.legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=2, fontsize=11)
plt.show()


# By UAS type and environment
fig, ax = plt.subplots(figsize=(6, 4.65))
vehicleCats = ["No UAS", "H520", "M300", "T150"]
data = dataByStimTestA.loc[:, ['AmbientEnv', 'UASType', 'NoticedPropFilt']]
data = data.sort_values('AmbientEnv', ignore_index=True)
data['UASType'] = pd.Categorical(values=data['UASType'],
                                 categories=vehicleCats)

sns.violinplot(data=data, y='NoticedPropFilt', split=True, hue='AmbientEnv',
               x='UASType', cut=0, palette=mycolours[0:2], inner='quart',
               width=0.65, bw_method='scott')

plt.setp(ax.collections, alpha=0.3)
# Add jittered dots
xjitter = 0.01
yjitter = 0.0
y_subdataPark = list()
y_subdataStreet = list()
for ii in range(0, len(vehicleCats)):
    y_subdataPark.append(dataByStimTestA.loc[np.logical_and(
                                             dataByStimTestA['AmbientEnv']
                                             == 'Park',
                                             dataByStimTestA['UASType']
                                             == vehicleCats[ii]),
                                             'NoticedPropFilt'].values)
    y_subdataStreet.append(dataByStimTestA.loc[np.logical_and(
                                               dataByStimTestA['AmbientEnv']
                                               == 'Street',
                                               dataByStimTestA['UASType']
                                               == vehicleCats[ii]),
                                               'NoticedPropFilt'].values)
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
    ax.scatter(x, y, s=10, color=mycolours[0], alpha=0.3)
for x, y in zip(x_Streetjittered, y_Streetjittered):
    ax.scatter(x, y, s=10, color=mycolours[1], alpha=0.3)

ax.set(yticks=np.arange(0, 1.1, 0.1), yticklabels=range(0, 110, 10),
       xticks=range(0, len(vehicleCats)), xticklabels=vehicleCats,
       xlabel="UAS type",
       ylabel="Proportion UAS noticed, %", ylim=[-0.05, 1.05])
ax.legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=2, fontsize=11)
plt.show()

# By flight operation and environment
fig, ax = plt.subplots(figsize=(6, 4.65))
opCats = ["No UAS", "Flyby", "Landing", "Takeoff"]
data = dataByStimTestA.loc[:, ['AmbientEnv', 'UASOperation', 'NoticedPropFilt']]
data = data.sort_values('AmbientEnv', ignore_index=True)
data['UASOperation'] = pd.Categorical(values=data['UASOperation'],
                                 categories=opCats)

sns.violinplot(data=data, y='NoticedPropFilt', split=True, hue='AmbientEnv',
               x='UASOperation', cut=0, palette=mycolours[0:2], inner='quart',
               width=0.65, bw_method='scott')

plt.setp(ax.collections, alpha=0.3)
# Add jittered dots
xjitter = 0.01
yjitter = 0.0
y_subdataPark = list()
y_subdataStreet = list()
for ii in range(0, len(opCats)):
    y_subdataPark.append(dataByStimTestA.loc[np.logical_and(
                                             dataByStimTestA['AmbientEnv']
                                             == 'Park',
                                             dataByStimTestA['UASOperation']
                                             == opCats[ii]),
                                             'NoticedPropFilt'].values)
    y_subdataStreet.append(dataByStimTestA.loc[np.logical_and(
                                               dataByStimTestA['AmbientEnv']
                                               == 'Street',
                                               dataByStimTestA['UASOperation']
                                               == opCats[ii]),
                                               'NoticedPropFilt'].values)
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
    ax.scatter(x, y, s=10, color=mycolours[0], alpha=0.3)
for x, y in zip(x_Streetjittered, y_Streetjittered):
    ax.scatter(x, y, s=10, color=mycolours[1], alpha=0.3)

ax.set(yticks=np.arange(0, 1.1, 0.1), yticklabels=range(0, 110, 10),
       xticks=range(0, len(opCats)), xticklabels=opCats,
       xlabel="UAS operation",
       ylabel="Proportion UAS noticed, %", ylim=[-0.05, 1.05])
ax.legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=2, fontsize=11)
plt.show()

# remove added columns from unfiltered dataset
dataByStimTestA.drop(labels=['NoticedPropFilt', 'NoticedTotalFilt'], axis=1,
                     inplace=True)




# -----------------------------------------------------------------------------
# Association analysis: Scatter plots & correlation
# -----------------------------------------------------------------------------

# reduce dataset for response analysis
includecols = ["LAeqMaxLR", "LAEMaxLR", "LAFmaxMaxLR",
               "LAF5ExMaxLR", "LAF10ExMaxLR", "LAF25ExMaxLR", "LAF50ExMaxLR",
               "LAF75ExMaxLR", "LAF90ExMaxLR", "LAF95ExMaxLR", "LASmaxMaxLR",
               "EPNLMaxLR", "EPNLNoTMaxLR", "PNLTmaxMaxLR",
               "LoudECMAHMSPowAvgBin", "TonalECMAHMSAvgMaxLR",
               "RoughECMAHMS10ExBin", "FluctstrHMS10ExBin",
               "LoudISO105ExMaxLR", "LoudISO1PowAvgMaxLR",
               "LoudISO3PowAvgMaxLR", "SharpAuresISO3AvgMaxLR",
               "SharpAuresISO3PowAvgMaxLR", "ImpulsHMSPowAvgMaxLR",
               "ImpulsHMSAvgMaxLR", "ImpulsHMSMaxMaxLR"]
includecols = (includecols +
               ["UAS" + includecol for includecol in includecols]
               + ["Amb" + includecol for includecol in includecols])
includecols = (["AnnoyMedian", "AnnoyMean", "ValenceMedian", "ValenceMean",
                "ArousalMedian", "ArousalMean"]
               + includecols)

ptADataForCorr = dataByStimTestA.loc[:, includecols]

UASTypes = dataByStimTestA['UASType'].unique()[1:]
ptADataForCorrByUAS = [dataTemp for dataTemp in
                       [dataByStimTestA.loc[dataByStimTestA['UASType'] == UAS,
                                            includecols] for UAS
                        in UASTypes]]

UASOps = dataByStimTestA['UASOperation'].unique()[1:]
ptADataForCorrByOp = [dataTemp for dataTemp in
                      [dataByStimTestA.loc[dataByStimTestA['UASOperation']
                                           == Op, includecols] for Op
                       in UASOps]]

ptADataForCorrByUASByOp = list()
count = 0
for UAS in UASTypes:
    for Op in UASOps:
        ptADataForCorrByUASByOp.append(dataByStimTestA.loc[np.logical_and(dataByStimTestA['UASType']
                                                                          == UAS,
                                                                          dataByStimTestA['UASOperation']
                                                                          == Op),
                                                           includecols])


# run correlation analysis
corr_matrixA = [ptADataForCorr.corr(method=corr_method, numeric_only=False)
                for corr_method in ['pearson', 'spearman', 'kendall']]
corr_matrixAByUAS = [[dataByUAS.corr(method=corr_method,
                                     numeric_only=False)
                      for corr_method in ['pearson', 'spearman', 'kendall']]
                     for dataByUAS in ptADataForCorrByUAS]
corr_matrixAByOp = [[dataByOp.corr(method=corr_method,
                                   numeric_only=False)
                     for corr_method in ['pearson', 'spearman', 'kendall']]
                    for dataByOp in ptADataForCorrByOp]
corr_matrixAByUASByOp = [[dataByUASByOp.corr(method=corr_method,
                                             numeric_only=False)
                          for corr_method in ['pearson', 'spearman',
                                              'kendall']]
                         for dataByUASByOp in ptADataForCorrByUASByOp]


# re-run correlation analysis to get p-values
corr_listA_pr = pd.DataFrame()
corr_listA_sp = pd.DataFrame()
corr_listA_kd = pd.DataFrame()
feat1s = [None]*len(ptADataForCorr.columns)**2
feat2s = [None]*len(ptADataForCorr.columns)**2
corrs_pr = np.zeros(len(ptADataForCorr.columns)**2)
corrs_sp = np.zeros(len(ptADataForCorr.columns)**2)
corrs_kd = np.zeros(len(ptADataForCorr.columns)**2)
pvalues_pr = np.zeros(len(ptADataForCorr.columns)**2)
pvalues_sp = np.zeros(len(ptADataForCorr.columns)**2)
pvalues_kd = np.zeros(len(ptADataForCorr.columns)**2)

count = 0
for feat1 in ptADataForCorr.columns:
    for feat2 in ptADataForCorr.columns:
        if feat1 == feat2:
            corr_pr, pvalue_pr = 1, 0
            corr_kd, pvalue_kd = 1, 0
        else:
            # avoid nan issues and speed up Kendall's tau
            mask = ~np.logical_or(ptADataForCorr[feat1].isna(),
                                  ptADataForCorr[feat2].isna())

            corr_pr, pvalue_pr = stats.pearsonr(ptADataForCorr.loc[mask,
                                                                   feat1].astype(float),
                                                ptADataForCorr.loc[mask,
                                                                   feat2].astype(float))
            corr_kd, pvalue_kd = stats.kendalltau(ptADataForCorr.loc[mask,
                                                                     feat1].astype(float),
                                                  ptADataForCorr.loc[mask,
                                                                     feat2].astype(float))
        corr_sp, pvalue_sp = stats.spearmanr(ptADataForCorr[feat1].astype(float),
                                             ptADataForCorr[feat2].astype(float),
                                             nan_policy='omit')

        feat1s[count] = feat1
        feat2s[count] = feat2
        corrs_pr[count] = corr_pr
        pvalues_pr[count] = pvalue_pr
        corrs_sp[count] = corr_sp
        pvalues_sp[count] = pvalue_sp
        corrs_kd[count] = corr_kd
        pvalues_kd[count] = pvalue_kd

        count += 1

# end of for loop for correlation analysis

corr_listA_pr['Variable_1'] = feat1s
corr_listA_pr['Variable_2'] = feat2s
corr_listA_pr['Correlation'] = corrs_pr
corr_listA_pr['p_value'] = pvalues_pr

corr_listA_sp['Variable_1'] = feat1s
corr_listA_sp['Variable_2'] = feat2s
corr_listA_sp['Correlation'] = corrs_sp
corr_listA_sp['p_value'] = pvalues_sp

corr_listA_kd['Variable_1'] = feat1s
corr_listA_kd['Variable_2'] = feat2s
corr_listA_kd['Correlation'] = corrs_kd
corr_listA_kd['p_value'] = pvalues_kd

corr_matrixA_p = [corrlist.pivot(columns='Variable_2',
                                 index='Variable_1',
                                 values='p_value')[corr_matrixA[1].columns].reindex(corr_matrixA[1].index)
                  for corrlist in [corr_listA_pr,
                                   corr_listA_sp,
                                   corr_listA_kd]]

# reduce dataset for noticeability  analysis
includecols = ["LAeqMaxLR", "LAEMaxLR", "LAFmaxMaxLR",
               "LAF5ExMaxLR", "LAF10ExMaxLR", "LAF25ExMaxLR", "LAF50ExMaxLR",
               "LAF75ExMaxLR", "LAF90ExMaxLR", "LAF95ExMaxLR", "LASmaxMaxLR",
               "EPNLMaxLR", "EPNLNoTMaxLR", "PNLTmaxMaxLR",
               "LoudECMAHMSPowAvgBin", "TonalECMAHMSAvgMaxLR",
               "RoughECMAHMS10ExBin", "FluctstrHMS10ExBin",
               "LoudISO105ExMaxLR", "LoudISO1PowAvgMaxLR",
               "LoudISO3PowAvgMaxLR", "SharpAuresISO3AvgMaxLR",
               "SharpAuresISO3PowAvgMaxLR", "ImpulsHMSPowAvgMaxLR",
               "ImpulsHMSAvgMaxLR", "ImpulsHMSMaxMaxLR"]
includecols = (includecols +
               ["UAS" + includecol for includecol in includecols]
               + ["Amb" + includecol for includecol in includecols])
includecols = includecols + ["SNRlevel", "UASLAeqdiffAmbLAF90",
                             "UASLASmaxdiffAmbLAF90", "UASLASmaxdiffAmbLAF50",
                             "UASLASmaxdiffAmbLAeq", "UASPartLoudGMSTPowAvg",
                             "UASPartLoudGMLTPowAvg",
                             "EPNLdiff", "EPNLNoTdiff", "PNLTmaxdiff"]
includecols = (["AnnoyMedianFilt", "AnnoyMeanFilt", "ValenceMedianFilt", "ValenceMeanFilt",
                "ArousalMedianFilt", "ArousalMeanFilt", "NoticedPropFilt"]
               + includecols)

ptANoticeDataForCorr = dataByStimTestANotice.loc[:, includecols]
# convert categorical data to numerical
ptANoticeDataForCorr['SNRlevel'] = dataByStimTestANotice['SNRlevel'].cat.codes
ptANoticeDataForCorr.loc[ptANoticeDataForCorr.loc[:, 'SNRlevel']
                         == -1, 'SNRlevel'] = np.nan
ptANoticeDataForCorr['SNRlevel'] = ptANoticeDataForCorr['SNRlevel'].astype(float)

UASTypes = dataByStimTestANotice['UASType'].unique()[1:]
ptANoticeDataForCorrByUAS = [dataTemp for dataTemp in
                             [dataByStimTestANotice.loc[dataByStimTestANotice['UASType'] == UAS,
                                                        includecols] for UAS
                              in UASTypes]]

UASOps = dataByStimTestANotice['UASOperation'].unique()[1:]
ptANoticeDataForCorrByOp = [dataTemp for dataTemp in
                            [dataByStimTestANotice.loc[dataByStimTestANotice['UASOperation']
                                                       == Op, includecols] for Op
                             in UASOps]]

ptANoticeDataForCorrByUASByOp = list()
count = 0
for UAS in UASTypes:
    for Op in UASOps:
        ptANoticeDataForCorrByUASByOp.append(dataByStimTestANotice.loc[np.logical_and(dataByStimTestANotice['UASType']
                                                                                 == UAS,
                                                                                 dataByStimTestANotice['UASOperation']
                                                                                 == Op),
                                                                       includecols])

# run correlation analysis
corr_matrixANotice = [ptANoticeDataForCorr.corr(method=corr_method, numeric_only=False)
                      for corr_method in ['pearson', 'spearman', 'kendall']]
corr_matrixANoticeByUAS = [[dataByUAS.corr(method=corr_method,
                                           numeric_only=False)
                            for corr_method in ['pearson', 'spearman', 'kendall']]
                           for dataByUAS in ptANoticeDataForCorrByUAS]
corr_matrixANoticeByOp = [[dataByOp.corr(method=corr_method,
                                   numeric_only=False)
                     for corr_method in ['pearson', 'spearman', 'kendall']]
                    for dataByOp in ptANoticeDataForCorrByOp]
corr_matrixANoticeByUASByOp = [[dataByUASByOp.corr(method=corr_method,
                                                   numeric_only=False)
                                for corr_method in ['pearson', 'spearman',
                                                    'kendall']]
                               for dataByUASByOp in ptANoticeDataForCorrByUASByOp]


# re-run correlation analysis to get p-values
corr_listA_pr = pd.DataFrame()
corr_listA_sp = pd.DataFrame()
corr_listA_kd = pd.DataFrame()
feat1s = [None]*len(ptANoticeDataForCorr.columns)**2
feat2s = [None]*len(ptANoticeDataForCorr.columns)**2
corrs_pr = np.zeros(len(ptANoticeDataForCorr.columns)**2)
corrs_sp = np.zeros(len(ptANoticeDataForCorr.columns)**2)
corrs_kd = np.zeros(len(ptANoticeDataForCorr.columns)**2)
pvalues_pr = np.zeros(len(ptANoticeDataForCorr.columns)**2)
pvalues_sp = np.zeros(len(ptANoticeDataForCorr.columns)**2)
pvalues_kd = np.zeros(len(ptANoticeDataForCorr.columns)**2)

count = 0
for feat1 in ptANoticeDataForCorr.columns:
    for feat2 in ptANoticeDataForCorr.columns:
        if feat1 == feat2:
            corr_pr, pvalue_pr = 1, 0
            corr_kd, pvalue_kd = 1, 0
        else:
            # avoid nan issues and speed up Kendall's tau
            mask = ~np.logical_or(ptANoticeDataForCorr[feat1].isna(),
                                  ptANoticeDataForCorr[feat2].isna())

            corr_pr, pvalue_pr = stats.pearsonr(ptANoticeDataForCorr.loc[mask,
                                                                   feat1].astype(float),
                                                ptANoticeDataForCorr.loc[mask,
                                                                   feat2].astype(float))
            corr_kd, pvalue_kd = stats.kendalltau(ptANoticeDataForCorr.loc[mask,
                                                                     feat1].astype(float),
                                                  ptANoticeDataForCorr.loc[mask,
                                                                     feat2].astype(float))
        corr_sp, pvalue_sp = stats.spearmanr(ptANoticeDataForCorr[feat1].astype(float),
                                             ptANoticeDataForCorr[feat2].astype(float),
                                             nan_policy='omit')

        feat1s[count] = feat1
        feat2s[count] = feat2
        corrs_pr[count] = corr_pr
        pvalues_pr[count] = pvalue_pr
        corrs_sp[count] = corr_sp
        pvalues_sp[count] = pvalue_sp
        corrs_kd[count] = corr_kd
        pvalues_kd[count] = pvalue_kd

        count += 1

# end of for loop for correlation analysis

corr_listA_pr['Variable_1'] = feat1s
corr_listA_pr['Variable_2'] = feat2s
corr_listA_pr['Correlation'] = corrs_pr
corr_listA_pr['p_value'] = pvalues_pr

corr_listA_sp['Variable_1'] = feat1s
corr_listA_sp['Variable_2'] = feat2s
corr_listA_sp['Correlation'] = corrs_sp
corr_listA_sp['p_value'] = pvalues_sp

corr_listA_kd['Variable_1'] = feat1s
corr_listA_kd['Variable_2'] = feat2s
corr_listA_kd['Correlation'] = corrs_kd
corr_listA_kd['p_value'] = pvalues_kd

corr_matrixANotice_p = [corrlist.pivot(columns='Variable_2',
                                       index='Variable_1',
                                       values='p_value')[corr_matrixANotice[1].columns].reindex(corr_matrixANotice[1].index)
                        for corrlist in [corr_listA_pr,
                                         corr_listA_sp,
                                         corr_listA_kd]]

# Noticeability correlation analysis
# ----------------------------------

# difference metrics
indicesDiff = ['SNRlevel', 'UASPartLoudGMLTPowAvg',
               'UASLAeqdiffAmbLAF90',
               'UASLASmaxdiffAmbLAF90', 'UASLASmaxdiffAmbLAF50',
               'UASLASmaxdiffAmbLAeq', 'EPNLdiff', 'EPNLNoTdiff',
               'PNLTmaxdiff']

labelsDiff = ["UAS $L_\\mathrm{Aeq}$ - ambient $L_\\mathrm{Aeq}$, dB",
              "LT Part. loudness (M-G), sones",
              "UAS $L_\\mathrm{Aeq}$ - ambient $L_\\mathrm{AF90}$, dB",
              "UAS $L_\\mathrm{ASmax}$ - ambient $L_\\mathrm{AF90}$, dB",
              "UAS $L_\\mathrm{ASmax}$ - ambient $L_\\mathrm{AF50}$, dB",
              "UAS $L_\\mathrm{ASmax}$ - ambient $L_\\mathrm{Aeq}$, dB",
              "UAS EPNL - ambient EPNL, dB",
              "UAS EPNL$_\\mathrm{no\\ no tone}$ - ambient EPNL$_\\mathrm{no\\ no tone}$, dB",
              "UAS PNLT$_\\mathrm{max}$ - ambient PNLT$_\\mathrm{max}$, dB"]


fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(9.5, 9.5))
for ii, (index, ax) in enumerate(zip(indicesDiff, axs.ravel())):
    if ii > len(indicesDiff):
        pass
    else:
        sns.scatterplot(data=ptANoticeDataForCorr, x=index, y='NoticedPropFilt',
                        ax=ax, color=mycolours[ii], alpha=0.4)
        ax.set(xticks=range(int((np.min(ptANoticeDataForCorr[index])-1)/5)*5,
                            round(np.ceil(np.max(ptANoticeDataForCorr[index]))/5)*5 + 5,
                            5), yticks=np.arange(0, 1.1, 0.1),
               xlabel=labelsDiff[ii],
               ylabel="Proportion UAS noticed, %")
        ax.grid(alpha=0.15, linestyle='--')
plt.show()


fig, ax = plt.subplots(figsize=(7.5, 4))
data = pd.DataFrame(data=[corr_matrixANotice[1].loc['NoticedPropFilt',
                                              indicesDiff],
                          corr_matrixANotice[2].loc['NoticedPropFilt',
                                              indicesDiff]],
                    index=["Spearman", "Kendall"])
sns.barplot(data=data.melt(var_name="Type",
                           ignore_index=False).reset_index(names="Method"),
            x='Type', y='value', hue='Method', width=0.6,
            palette=mycolours[0:2], ax=ax)
labelsDiffShort = ["$\\Delta L_\\mathrm{Aeq}$",
                   "LT Part. loudness (M-G)",
                   "$\\Delta(L_\\mathrm{Aeq}-L_\\mathrm{AF90})$",
                   "$\\Delta(L_\\mathrm{ASmax}-L_\\mathrm{AF90})$",
                   "$\\Delta(L_\\mathrm{ASmax}-L_\\mathrm{AF50})$",
                   "$\\Delta(L_\\mathrm{ASmax}-L_\\mathrm{Aeq})$",
                   "$\\Delta($EPNL$)$",
                   "$\\Delta($EPNL$_\\mathrm{no\\ tone})$",
                   "$\\Delta($PNLT$_\\mathrm{max})$"]
ax.set(yticks=np.arange(0, 1.1, 0.1), xlabel="UAS sound vs ambient sound",
       ylabel="Correlation coefficient" + "\n" + "(proportion UAS noticed, %)")
ax.set_xticks(ticks=range(len(labelsDiffShort)), labels=labelsDiffShort,
              rotation=45, ha='right', rotation_mode='anchor')
ax.legend(loc='best')
ax.grid(alpha=0.15, linestyle='--', axis='y')
plt.show()


# By operation
UASOps = dataByStimTestA['UASOperation'].unique()[1:]
data = [pd.DataFrame(data=[corr_matrixANoticeByOp[0][1].loc['NoticedPropFilt',
                                                            indicesDiff],
                           corr_matrixANoticeByOp[1][1].loc['NoticedPropFilt',
                                                            indicesDiff],
                           corr_matrixANoticeByOp[2][1].loc['NoticedPropFilt',
                                                            indicesDiff]],
                     index=["Flyover", "Landing", "Takeoff"]),
        pd.DataFrame(data=[corr_matrixANoticeByOp[0][2].loc['NoticedPropFilt',
                                                            indicesDiff],
                           corr_matrixANoticeByOp[1][2].loc['NoticedPropFilt',
                                                            indicesDiff],
                           corr_matrixANoticeByOp[2][2].loc['NoticedPropFilt',
                                                            indicesDiff]],
                     index=["Flyover", "Landing", "Takeoff"])]

# fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7.5, 6))
# for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
#                                           axs.ravel())):
#     sns.heatmap(data=data[iiCoeff], vmin=0.35, vmax=0.85,
#                 xticklabels=labelsDiffShort, yticklabels=data[iiCoeff].index,
#                 ax=ax, cbar_kws=dict(label="Correlation coefficient"))
#     ax.set_xticks(ticks=range(len(labelsDiffShort)), labels=labelsDiffShort,
#                   rotation=30)
#     ax.set_title("  " + Coeff, y=1.0, pad=5, loc='left')
# plt.show()

fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(9, 6))
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.barplot(data=data[iiCoeff].melt(var_name="Type",
                                        ignore_index=False).reset_index(names="Operation"),
                x='Type', y='value', hue='Operation', width=0.5,
                palette=mycolours[0:3], ax=ax)
    if iiCoeff < len(["Spearman", "Kendall"]) - 1:
        ax.set_xticks(ticks=range(len(labelsDiffShort)), labels=[])
        ax.set(xlabel=None)
    else:
        ax.set(xlabel="UAS sound vs ambient sound")
        ax.set_xticks(ticks=range(len(labelsDiffShort)),
                      labels=labelsDiffShort, rotation=45,
                      ha='right', rotation_mode='anchor')

    ax.set(yticks=np.arange(0, 1.1, 0.1),
           ylabel="Correlation coefficient" + "\n" + "(proportion UAS noticed, %)")
    ax.set_title("  " + Coeff, y=1.0, pad=-14, loc='left')

    ax.legend(loc='upper center', ncol=len(UASOps))
    ax.grid(alpha=0.15, linestyle='--', axis='y')
plt.show()

# By UAS type
UASTypes = dataByStimTestA['UASType'].unique()[1:]
data = [pd.DataFrame(data=[corr_matrixANoticeByUAS[0][1].loc['NoticedPropFilt',
                                                       indicesDiff],
                           corr_matrixANoticeByUAS[1][1].loc['NoticedPropFilt',
                                                       indicesDiff],
                           corr_matrixANoticeByUAS[2][1].loc['NoticedPropFilt',
                                                       indicesDiff]],
                     index=["H520E", "M300", "T150"]),
        pd.DataFrame(data=[corr_matrixANoticeByUAS[0][2].loc['NoticedPropFilt',
                                                       indicesDiff],
                           corr_matrixANoticeByUAS[1][2].loc['NoticedPropFilt',
                                                       indicesDiff],
                           corr_matrixANoticeByUAS[2][2].loc['NoticedPropFilt',
                                                       indicesDiff]],
                     index=["H520E", "M300", "T150"])]

# fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7.5, 6))
# for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
#                                           axs.ravel())):
#     sns.heatmap(data=data[iiCoeff], vmin=0.35, vmax=0.85,
#                 xticklabels=labelsDiffShort, yticklabels=data[iiCoeff].index,
#                 ax=ax, cbar_kws=dict(label="Correlation coefficient"
#                                      + "\n" + "(UAS noticed, %)"))
#     ax.set_xticks(ticks=range(len(labelsDiffShort)), labels=labelsDiffShort,
#                   rotation=30)
#     ax.set_title("  " + Coeff, y=1.0, pad=5, loc='left')
# plt.show()

fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(9, 6))
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.barplot(data=data[iiCoeff].melt(var_name="Type",
                                        ignore_index=False).reset_index(names="UAS Type"),
                x='Type', y='value', hue='UAS Type', width=0.5,
                palette=mycolours[0:3], ax=ax)
    if iiCoeff < len(["Spearman", "Kendall"]) - 1:
        ax.set_xticks(ticks=range(len(labelsDiffShort)), labels=[])
        ax.set(xlabel=None)
    else:
        ax.set(xlabel="UAS sound vs ambient sound")
        ax.set_xticks(ticks=range(len(labelsDiffShort)),
                      labels=labelsDiffShort, rotation=45,
                      ha='right', rotation_mode='anchor')

    ax.set(yticks=np.arange(0, 1.1, 0.1),
           ylabel="Correlation coefficient" + "\n" + "(proportion UAS noticed, %)")
    ax.set_title("  " + Coeff, y=1.0, pad=-14, loc='left')

    ax.legend(loc='upper center', ncol=len(UASTypes))
    ax.grid(alpha=0.15, linestyle='--', axis='y')
plt.show()

# sound quality metrics
indicesSonQual = ['TonalECMAHMSAvgMaxLR', 'RoughECMAHMS10ExBin',
                  'FluctstrHMS10ExBin', 'SharpAuresISO3AvgMaxLR',
                  'SharpAuresISO3PowAvgMaxLR', 'ImpulsHMSAvgMaxLR',
                  'ImpulsHMSMaxMaxLR', 'ImpulsHMSPowAvgMaxLR']
indicesSonQual = indicesSonQual + (["UAS" + index
                                    for index in indicesSonQual]
                                   + ["Amb" + index
                                      for index in indicesSonQual])

labelsSonQual = ["$T_\\mathrm{mean,ECMA(HMS)}$", "$R_\\mathrm{10,ECMA(HMS)}$",
                 "$F_\\mathrm{10,HMS}$", "$S_\\mathrm{mean,Aures,ISO3}$",
                 "$S_\\mathrm{Aures,ISO3}$", "$I_\\mathrm{mean,HMS}$",
                 "$I_\\mathrm{max,HMS}$", "$I_\\mathrm{HMS}$"]
labelsSonQual = labelsSonQual + (["UAS " + label for label in labelsSonQual]
                                 + ["Amb " + label for label in labelsSonQual])

fig, ax = plt.subplots(figsize=(8, 4))
data = pd.DataFrame(data=[corr_matrixANotice[1].loc['NoticedPropFilt', indicesSonQual],
                          corr_matrixANotice[2].loc['NoticedPropFilt', indicesSonQual]],
                    index=["Spearman", "Kendall"])
sns.barplot(data=data.melt(var_name="Type",
                           ignore_index=False).reset_index(names="Method"),
            x='Type', y='value', hue='Method', width=0.75,
            palette=mycolours[0:2], ax=ax)
ax.set(yticks=np.arange(-0.6, 1.1, 0.1), xlabel="Sound quality metric",
       ylabel="Correlation coefficient" + "\n" + "(proportion UAS noticed)")
ax.set_xticks(ticks=range(len(labelsSonQual)), labels=labelsSonQual,
              rotation=45, ha='right', rotation_mode='anchor')
ax.legend(loc='best')
ax.grid(alpha=0.15, linestyle='--', axis='y')
plt.show()

# plot UAS sound quality data
labelsSonQualSub = (labelsSonQual[8:11] + [labelsSonQual[12]]
                    + [labelsSonQual[15]])
fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(15, 3.75))
for ii, (index, ax) in enumerate(zip(indicesSonQual[8:11]
                                     + [indicesSonQual[12]]
                                     + [indicesSonQual[15]],
                                     axs.ravel())):
    sns.regplot(data=ptANoticeDataForCorr, x=index, y='AnnoyMedianFilt',
                ax=ax, y_jitter=0.1, color=mycolours[ii],
                scatter=True, fit_reg=False, scatter_kws=dict(alpha=0.4))
    if index == 'UASFluctstrHMS10ExBin':
        ax.set(xticks=np.arange(np.round(np.min(ptANoticeDataForCorr[index]), 1),
                                np.ceil(np.max(ptANoticeDataForCorr[index])*100)/100,
                                np.ceil(np.max(ptANoticeDataForCorr[index])*100)/100/5))
    else:
        ax.set(xticks=np.arange(np.round(np.min(ptANoticeDataForCorr[index]), 2),
                                np.round(np.ceil(np.max(ptANoticeDataForCorr[index])*10)/10, 1),
                                np.round(np.ceil(np.max(ptANoticeDataForCorr[index])*10)/10, 1)/5))
    ax.set(yticks=range(0, 11),
           xlabel=labelsSonQualSub[ii],
           ylabel="Median annoyance rating")
    ax.grid(alpha=0.15, linestyle='--')
plt.show()

# plot whole scene sound quality data
labelsSonQualSub = (labelsSonQual[0:3] + [labelsSonQual[4]]
                    + [labelsSonQual[6]])
fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(15, 3.75))
ptADataForCorrSub1 = ptADataForCorr.loc[~ptADataForCorr.index.isin(["A1_CALBIN_Pa.wav",
                                                                    "A2_CALBIN_Pa.wav"])]
ptADataForCorrSub2 = ptADataForCorr.loc[ptADataForCorr.index.isin(["A1_CALBIN_Pa.wav",
                                                                   "A2_CALBIN_Pa.wav"])]
for ii, (index, ax) in enumerate(zip(indicesSonQual[0:3]
                                     + [indicesSonQual[4]]
                                     + [indicesSonQual[6]],
                                     axs.ravel())):
    sns.regplot(data=ptADataForCorrSub1, x=index, y='AnnoyMedian',
                ax=ax, y_jitter=0.1, color=mycolours[ii],
                scatter=True, fit_reg=False, scatter_kws=dict(alpha=0.4))
    sns.regplot(data=ptADataForCorrSub2, x=index, y='AnnoyMedian',
                ax=ax, y_jitter=0.1, color=mycolours[ii], marker='x',
                scatter=True, fit_reg=False, scatter_kws=dict(alpha=0.4))
    if index == 'UASFluctstrHMS10ExBin':
        ax.set(xticks=np.arange(np.round(np.min(ptADataForCorr[index]), 1),
                                np.ceil(np.max(ptADataForCorr[index])*100)/100,
                                np.ceil(np.max(ptADataForCorr[index])*100)/100/5))
    else:
        ax.set(xticks=np.arange(np.round(np.min(ptADataForCorr[index]), 2),
                                np.round(np.ceil(np.max(ptADataForCorr[index])*10)/10, 1),
                                np.round(np.ceil(np.max(ptADataForCorr[index])*10)/10, 1)/5))
    ax.set(yticks=range(0, 11),
           xlabel=labelsSonQualSub[ii],
           ylabel="Median annoyance rating")
    ax.grid(alpha=0.15, linestyle='--')
plt.show()


# By operation
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(13, 6))
data = [pd.DataFrame(data=[corr_matrixAByOp[0][1].loc['AnnoyMedian',
                                                      indicesSonQual],
                           corr_matrixAByOp[1][1].loc['AnnoyMedian',
                                                      indicesSonQual],
                           corr_matrixAByOp[2][1].loc['AnnoyMedian',
                                                      indicesSonQual]],
                     index=["Flyover", "Landing", "Takeoff"]),
        pd.DataFrame(data=[corr_matrixAByOp[0][2].loc['AnnoyMedian',
                                                      indicesSonQual],
                           corr_matrixAByOp[1][2].loc['AnnoyMedian',
                                                      indicesSonQual],
                           corr_matrixAByOp[2][2].loc['AnnoyMedian',
                                                      indicesSonQual]],
                     index=["Flyover", "Landing", "Takeoff"])]
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.heatmap(data=data[iiCoeff], vmin=-0.5, vmax=0.7,
                xticklabels=labelsSonQual, yticklabels=data[iiCoeff].index,
                ax=ax, cmap='seismic',
                cbar_kws=dict(label="Correlation coefficient"
                              + "\n" + "(median annoyance rating)"))
    ax.set_xticks(ticks=range(len(labelsSonQual)), labels=labelsSonQual,
                  rotation=30)
    ax.set_title("  " + Coeff, y=1.0, pad=5, loc='left')


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(13, 6))
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.barplot(data=data[iiCoeff].melt(var_name="Type",
                                        ignore_index=False).reset_index(names="Operation"),
                x='Type', y='value', hue='Operation', width=0.5,
                palette=mycolours[0:3], ax=ax)
    if iiCoeff < len(["Spearman", "Kendall"]) - 1:
        ax.set_xticks(ticks=range(len(labelsSonQual)), labels=[])
        ax.set(xlabel=None)
    else:
        ax.set(xlabel="Sound quality metric")
        ax.set_xticks(ticks=range(len(labelsSonQual)), labels=labelsSonQual,
                      rotation=45, ha='right', rotation_mode='anchor')

    ax.set(yticks=np.arange(-0.6, 1.1, 0.1),
           ylabel="Correlation coefficient"
           + "\n" + "(median annoyance rating)")
    ax.set_title("  " + Coeff, y=1.0, pad=-10, loc='left')

    ax.legend(loc='upper right', ncol=len(UASOps))
    ax.grid(alpha=0.15, linestyle='--', axis='y')


# By UAS type
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(13, 6))
data = [pd.DataFrame(data=[corr_matrixAByUAS[0][1].loc['AnnoyMedian',
                                                       indicesSonQual],
                           corr_matrixAByUAS[1][1].loc['AnnoyMedian',
                                                       indicesSonQual],
                           corr_matrixAByUAS[2][1].loc['AnnoyMedian',
                                                       indicesSonQual]],
                     index=["H520E", "M300", "T150"]),
        pd.DataFrame(data=[corr_matrixAByUAS[0][2].loc['AnnoyMedian',
                                                       indicesSonQual],
                           corr_matrixAByUAS[1][2].loc['AnnoyMedian',
                                                       indicesSonQual],
                           corr_matrixAByUAS[2][2].loc['AnnoyMedian',
                                                       indicesSonQual]],
                     index=["H520E", "M300", "T150"])]
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.heatmap(data=data[iiCoeff], vmin=-0.5, vmax=0.8,
                xticklabels=labelsSonQual, yticklabels=data[iiCoeff].index,
                ax=ax, cmap='seismic',
                cbar_kws=dict(label="Correlation coefficient"
                              + "\n" + "(median annoyance rating)"))
    ax.set_xticks(ticks=range(len(labelsSonQual)), labels=labelsSonQual,
                  rotation=30)
    ax.set_title("  " + Coeff, y=1.0, pad=5, loc='left')


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(13, 6))
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.barplot(data=data[iiCoeff].melt(var_name="Type",
                                        ignore_index=False).reset_index(names="UAS Type"),
                x='Type', y='value', hue='UAS Type', width=0.5,
                palette=mycolours[0:3], ax=ax)
    if iiCoeff < len(["Spearman", "Kendall"]) - 1:
        ax.set_xticks(ticks=range(len(labelsSonQual)), labels=[])
        ax.set(xlabel=None)
    else:
        ax.set(xlabel="Sound quality metric")
        ax.set_xticks(ticks=range(len(labelsSonQual)), labels=labelsSonQual,
                      rotation=45, ha='right', rotation_mode='anchor')

    ax.set(yticks=np.arange(-0.6, 1.1, 0.1),
           ylabel="Correlation coefficient"
           + "\n" + "(median annoyance rating)")
    ax.set_title("  " + Coeff, y=1.0, pad=-10, loc='left')

    ax.legend(loc='upper right', ncol=len(UASTypes))
    ax.grid(alpha=0.15, linestyle='--', axis='y')


# correlation notice vs median annoyance
fig, ax = plt.subplots(figsize=(5, 4))
data = dataByStimTestA.loc[:, ['AmbientEnv', 'NoticedProportion',
                               'AnnoyMedian']].sort_values(by='AmbientEnv')
sns.regplot(data=data.loc[data['AmbientEnv'] == "Park"],
            x='NoticedProportion', y='AnnoyMedian',
            color=mycolours[0], fit_reg=False, scatter=True, marker='o',
            scatter_kws=dict(alpha=0.4), y_jitter=0.1, ax=ax)
sns.regplot(data=data.loc[data['AmbientEnv'] == "Street"],
            x='NoticedProportion', y='AnnoyMedian',
            color=mycolours[1], fit_reg=False, scatter=True, marker='^',
            scatter_kws=dict(alpha=0.4), y_jitter=0.1, ax=ax)
ax.set(xticks=np.arange(0, 1.1, 0.1), xticklabels=range(0, 110, 10),
       xlim=[-0.05, 1.05], xlabel="Proportion UAS noticed, %",
       yticks=range(0, 11), ylabel="Median annoyance rating")
ax.grid(alpha=0.15, linestyle='--')
ax.legend(labels=["Park", "Street"], title="Ambient environment", loc='best')
plt.show()

# fig, ax = plt.subplots(figsize=(9, 4))
# data = dataByStimTestA.loc[:, ['AmbientEnv', 'UASOperation',
#                                'UASType', 'NoticedProportion',
#                                'AnnoyMedian']].sort_values(by='AmbientEnv')
# UASTypes = dataByStimTestA['UASType'].unique()[1:]
# UASOps = dataByStimTestA['UASOperation'].unique()[1:]
# markShapes = ['o', 's', '^']
# markEdge = mycolours[2:]
# labels = list()

# for ii, AmbEnv in enumerate(data.AmbientEnv.unique()):
#     # plot No UAS points
#     sns.regplot(data=data.loc[np.logical_and(data['AmbientEnv']
#                                              == AmbEnv,
#                                              data['UASType']
#                                              == "No UAS",
#                                              data['UASOperation']
#                                              == "No UAS")],
#                 x='NoticedProportion', y='AnnoyMedian',
#                 scatter=True, fit_reg=False, color=mycolours[ii],
#                 marker='+', y_jitter=0.1,
#                 scatter_kws=dict(alpha=0.3))
#     labels.append(AmbEnv + ": No UAS")
#     # plot remaining points
#     for jj, UASType in enumerate(UASTypes):
#         for kk, UASOp in enumerate(UASOps):
#             sns.regplot(data=data.loc[np.logical_and(data['AmbientEnv']
#                                                      == AmbEnv,
#                                                      data['UASType']
#                                                      == UASType,
#                                                      data['UASOperation']
#                                                      == UASOp)],
#                         x='NoticedProportion', y='AnnoyMedian',
#                         scatter=True, fit_reg=False, color=mycolours[ii],
#                         marker=markShapes[jj], y_jitter=0.1,
#                         scatter_kws=dict(edgecolors=markEdge[kk], alpha=0.3))
#             labels.append(AmbEnv + ": " + UASType + ", " + UASOp)
# ax.set(xticks=np.arange(0, 1.1, 0.1), xticklabels=range(0, 110, 10),
#        xlim=[-0.05, 1.05], xlabel="Proportion UAS noticed, %",
#        yticks=range(0, 11), ylabel="Median annoyance rating")
# ax.legend(labels=labels, ncol=2, bbox_to_anchor=(1.05, 1.0), loc='upper left')
# ax.grid(alpha=0.15, linestyle='--')
# plt.show()

fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(12, 9))
data = dataByStimTestA.loc[:, ['AmbientEnv', 'UASOperation',
                               'UASType', 'NoticedProportion',
                               'AnnoyMedian']].sort_values(by='AmbientEnv')
markShapes = ['o', '^']
dataForLoop = list()
for UASType in UASTypes:
    for UASOp in UASOps:
        dataForLoop.append(data.loc[np.logical_and(data['UASType']
                                                   == UASType,
                                                   data['UASOperation']
                                                   == UASOp)])

for ii, ax in enumerate(axs.ravel()):
    for jj, AmbEnv in enumerate(data.AmbientEnv.unique()):
        sns.regplot(data=data.loc[np.logical_and(data['AmbientEnv']
                                                 == AmbEnv,
                                                 data['UASType']
                                                 == "No UAS",
                                                 data['UASOperation']
                                                 == "No UAS")],
                    x='NoticedProportion', y='AnnoyMedian',
                    scatter=True, fit_reg=False, color=mycolours[jj],
                    marker='x', y_jitter=0.1,
                    scatter_kws=dict(alpha=0.4), ax=ax)
        sns.regplot(data=dataForLoop[ii].loc[data['AmbientEnv'] == AmbEnv],
                    x='NoticedProportion', y='AnnoyMedian',
                    scatter=True, fit_reg=False, color=mycolours[jj],
                    y_jitter=0.1, marker=markShapes[jj],
                    scatter_kws=dict(alpha=0.4), ax=ax)
    ax.set(xticks=np.arange(0, 1.1, 0.1), xticklabels=range(0, 110, 10),
           xlim=[-0.05, 1.05], xlabel="Proportion UAS noticed, %",
           yticks=range(0, 11), ylabel="Median annoyance rating")
    ax.grid(alpha=0.15, linestyle='--')
    ax.set_title(dataForLoop[ii].loc[:, 'UASType'].iloc[0]
                 + ": " + dataForLoop[ii].loc[:, 'UASOperation'].iloc[0],
                 y=1, pad=-15, loc='center')

plt.show()


# response correlation analysis
# -----------------------------

# plot annoyance correlations
# sound level indices
indicesAcoustic = ['LAeqMaxLR', 'LAFmaxMaxLR', 'LAF5ExMaxLR', 'LAF10ExMaxLR',
                   'LAF25ExMaxLR', 'LAF50ExMaxLR', 'LAF75ExMaxLR',
                   'LAF90ExMaxLR', 'LAF95ExMaxLR', 'LASmaxMaxLR', 'EPNLMaxLR',
                   'EPNLNoTMaxLR', 'PNLTmaxMaxLR']
indicesAcoustic = indicesAcoustic + (["UAS" + index
                                      for index in indicesAcoustic]
                                     + ["Amb" + index
                                        for index in indicesAcoustic])

labelsAcoustic = ["$L_\\mathrm{Aeq}$", "$L_\\mathrm{AFmax}$",
                  "$L_\\mathrm{AF5}$", "$L_\\mathrm{AF10}$",
                  "$L_\\mathrm{AF25}$", "$L_\\mathrm{AF50}$",
                  "$L_\\mathrm{AF75}$", "$L_\\mathrm{AF90}$",
                  "$L_\\mathrm{AF95}$", "$L_\\mathrm{ASmax}$",
                  "EPNL", "EPNL,$_\\mathrm{no\\ tone}$",
                  "PNLT$_\\mathrm{max}$"]
labelsAcoustic = labelsAcoustic + (["UAS " + label
                                    for label in labelsAcoustic]
                                   + ["Amb " + label
                                      for label in labelsAcoustic])

data = pd.DataFrame(data=[corr_matrixA[1].loc['AnnoyMedian', indicesAcoustic],
                          corr_matrixA[2].loc['AnnoyMedian', indicesAcoustic]],
                    index=["Spearman", "Kendall"])

fig, ax = plt.subplots(figsize=(13, 4))
sns.barplot(data=data.melt(var_name="Type",
                           ignore_index=False).reset_index(names="Method"),
            x='Type', y='value', hue='Method', width=0.75,
            palette=mycolours[0:2], ax=ax)
ax.set(yticks=np.arange(-0.5, 1.1, 0.1), xlabel="Sound level index",
       ylabel="Correlation coefficient" + "\n" + "(median annoyance rating)")
ax.set_xticks(ticks=range(len(labelsAcoustic)), labels=labelsAcoustic,
              rotation=45, ha='right', rotation_mode='anchor')
ax.legend(loc='best', ncol=2)
ax.grid(alpha=0.15, linestyle='--', axis='y')
plt.show()

fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(14, 7))
labelsAcousticSub = labelsAcoustic[13:17] + labelsAcoustic[22:26]
for ii, (index, ax) in enumerate(zip(indicesAcoustic[13:17]
                                     + indicesAcoustic[22:26],
                                     axs.ravel())):
    sns.regplot(data=ptADataForCorr, x=index, y='AnnoyMedian',
                ax=ax, y_jitter=0.1, color=mycolours[ii],
                scatter=True, fit_reg=False, scatter_kws=dict(alpha=0.4))
    ax.set(xticks=range(int(np.min(ptADataForCorr[index])/5)*5,
                        round(np.ceil(np.max(ptADataForCorr[index]))/5)*5 + 5,
                        5), yticks=range(0, 11),
           xlabel=labelsAcousticSub[ii] + ", dB",
           ylabel="Median annoyance rating")
    ax.grid(alpha=0.15, linestyle='--')

# By operation
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(15, 6))
data = [pd.DataFrame(data=[corr_matrixAByOp[0][1].loc['AnnoyMedian',
                                                      indicesAcoustic],
                           corr_matrixAByOp[1][1].loc['AnnoyMedian',
                                                      indicesAcoustic],
                           corr_matrixAByOp[2][1].loc['AnnoyMedian',
                                                      indicesAcoustic]],
                     index=["Flyover", "Landing", "Takeoff"]),
        pd.DataFrame(data=[corr_matrixAByOp[0][2].loc['AnnoyMedian',
                                                      indicesAcoustic],
                           corr_matrixAByOp[1][2].loc['AnnoyMedian',
                                                      indicesAcoustic],
                           corr_matrixAByOp[2][2].loc['AnnoyMedian',
                                                      indicesAcoustic]],
                     index=["Flyover", "Landing", "Takeoff"])]
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.heatmap(data=data[iiCoeff], vmin=-0.4, vmax=0.9,
                xticklabels=labelsAcoustic, yticklabels=data[iiCoeff].index,
                ax=ax, cbar_kws=dict(label="Correlation coefficient"
                                     + "\n" + "(median annoyance rating)"),
                cmap='seismic')
    ax.set_xticks(ticks=range(len(labelsAcoustic)), labels=labelsAcoustic,
                  rotation=30, ha='right', rotation_mode='anchor')
    ax.set_title("  " + Coeff, y=1.0, pad=5, loc='left')


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(18, 6))
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.barplot(data=data[iiCoeff].melt(var_name="Type",
                                        ignore_index=False).reset_index(names="Operation"),
                x='Type', y='value', hue='Operation', width=0.5,
                palette=mycolours[0:3], ax=ax)
    if iiCoeff < len(["Spearman", "Kendall"]) - 1:
        ax.set_xticks(ticks=range(len(labelsAcoustic)), labels=[])
        ax.set(xlabel=None)
    else:
        ax.set(xlabel="Sound level index")
        ax.set_xticks(ticks=range(len(labelsAcoustic)), labels=labelsAcoustic,
                      rotation=45, ha='right', rotation_mode='anchor')

    ax.set(yticks=np.arange(-0.5, 1.1, 0.1),
           ylabel="Correlation coefficient" + "\n" + "(median annoyance rating)")
    ax.set_title("  " + Coeff, y=0, pad=5, loc='left')

    ax.legend(loc='upper right', ncol=len(UASOps))
    ax.grid(alpha=0.15, linestyle='--', axis='y')


# By UAS type
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(15, 6))
data = [pd.DataFrame(data=[corr_matrixAByUAS[0][1].loc['AnnoyMedian',
                                                       indicesAcoustic],
                           corr_matrixAByUAS[1][1].loc['AnnoyMedian',
                                                       indicesAcoustic],
                           corr_matrixAByUAS[2][1].loc['AnnoyMedian',
                                                       indicesAcoustic]],
                     index=["H520E", "M300", "T150"]),
        pd.DataFrame(data=[corr_matrixAByUAS[0][2].loc['AnnoyMedian',
                                                       indicesAcoustic],
                           corr_matrixAByUAS[1][2].loc['AnnoyMedian',
                                                       indicesAcoustic],
                           corr_matrixAByUAS[2][2].loc['AnnoyMedian',
                                                       indicesAcoustic]],
                     index=["H520E", "M300", "T150"])]
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.heatmap(data=data[iiCoeff], vmin=-0.4, vmax=0.9,
                xticklabels=labelsAcoustic, yticklabels=data[iiCoeff].index,
                ax=ax, cmap='seismic',
                cbar_kws=dict(label="Correlation coefficient"
                              + "\n" + "(median annoyance rating)"))
    ax.set_xticks(ticks=range(len(labelsAcoustic)), labels=labelsAcoustic,
                  rotation=30, ha='right', rotation_mode='anchor')
    ax.set_title("  " + Coeff, y=1.0, pad=5, loc='left')


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(18, 6))
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.barplot(data=data[iiCoeff].melt(var_name="Type",
                                        ignore_index=False).reset_index(names="UAS Type"),
                x='Type', y='value', hue='UAS Type', width=0.5,
                palette=mycolours[0:3], ax=ax)
    if iiCoeff < len(["Spearman", "Kendall"]) - 1:
        ax.set_xticks(ticks=range(len(labelsAcoustic)), labels=[])
        ax.set(xlabel=None)
    else:
        ax.set(xlabel="Sound level index")
        ax.set_xticks(ticks=range(len(labelsAcoustic)), labels=labelsAcoustic,
                      rotation=45, ha='right', rotation_mode='anchor')

    ax.set(yticks=np.arange(-0.5, 1.1, 0.1),
           ylabel="Correlation coefficient"
           + "\n" + "(median annoyance rating)")
    ax.set_title("  " + Coeff, y=0, pad=5, loc='left')

    ax.legend(loc='upper right', ncol=len(UASTypes))
    ax.grid(alpha=0.15, linestyle='--', axis='y')


# loudness metrics
indicesLoudness = ['LoudECMAHMSPowAvgBin', 'LoudISO105ExMaxLR',
                   'LoudISO1PowAvgMaxLR', 'LoudISO3PowAvgMaxLR']
indicesLoudness = indicesLoudness + (["UAS" + index
                                      for index in indicesLoudness]
                                     + ["Amb" + index
                                        for index in indicesLoudness])

labelsLoud = ["$N_\\mathrm{ECMA(HMS)}$", "$N_{5\\mathrm{,ISO1}}$",
              "$N_\\mathrm{ISO1}$", "$N_\\mathrm{ISO3}$"]
labelsLoud = labelsLoud + (["UAS " + label for label in labelsLoud]
                           + ["Amb " + label for label in labelsLoud])

fig, ax = plt.subplots(figsize=(5, 4))
data = pd.DataFrame(data=[corr_matrixA[1].loc['AnnoyMedian', indicesLoudness],
                          corr_matrixA[2].loc['AnnoyMedian', indicesLoudness]],
                    index=["Spearman", "Kendall"])
sns.barplot(data=data.melt(var_name="Type",
                           ignore_index=False).reset_index(names="Method"),
            x='Type', y='value', hue='Method', width=0.75,
            palette=mycolours[0:2], ax=ax)
ax.set(yticks=np.arange(-0.3, 1.1, 0.1), xlabel="Loudness metric",
       ylabel="Correlation coefficient" + "\n" + "(median annoyance rating)")
ax.set_xticks(ticks=range(len(labelsLoud)), labels=labelsLoud, rotation=45,
              ha='right', rotation_mode='anchor')
ax.legend(loc='best')
ax.grid(alpha=0.15, linestyle='--', axis='y')


labelsLoudSub = [labelsLoud[4]] + [labelsLoud[7]]
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(7, 3.75))
for ii, (index, ax) in enumerate(zip([indicesLoudness[4]]
                                     + [indicesLoudness[7]],
                                     axs.ravel())):
    sns.regplot(data=ptADataForCorr, x=index, y='AnnoyMedian',
                ax=ax, y_jitter=0.1, color=mycolours[ii],
                scatter=True, fit_reg=False, scatter_kws=dict(alpha=0.4))
    ax.set(xticks=range(int(np.min(ptADataForCorr[index])/2)*2,
                        round(np.ceil(np.max(ptADataForCorr[index]))/2)*2 + 2,
                        2), yticks=range(0, 11),
           xlabel=labelsLoudSub[ii] + ", sones",
           ylabel="Median annoyance rating")
    ax.grid(alpha=0.15, linestyle='--')

# By operation
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7.5, 6))
data = [pd.DataFrame(data=[corr_matrixAByOp[0][1].loc['AnnoyMedian',
                                                      indicesLoudness],
                           corr_matrixAByOp[1][1].loc['AnnoyMedian',
                                                      indicesLoudness],
                           corr_matrixAByOp[2][1].loc['AnnoyMedian',
                                                      indicesLoudness]],
                     index=["Flyover", "Landing", "Takeoff"]),
        pd.DataFrame(data=[corr_matrixAByOp[0][2].loc['AnnoyMedian',
                                                      indicesLoudness],
                           corr_matrixAByOp[1][2].loc['AnnoyMedian',
                                                      indicesLoudness],
                           corr_matrixAByOp[2][2].loc['AnnoyMedian',
                                                      indicesLoudness]],
                     index=["Flyover", "Landing", "Takeoff"])]
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.heatmap(data=data[iiCoeff], vmin=0.3, vmax=0.9,
                xticklabels=labelsLoud, yticklabels=data[iiCoeff].index,
                ax=ax, cbar_kws=dict(label="Correlation coefficient"
                                     + "\n" + "(median annoyance rating)"))
    ax.set_xticks(ticks=range(len(labelsLoud)), labels=labelsLoud,
                  rotation=30)
    ax.set_title("  " + Coeff, y=1.0, pad=5, loc='left')


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(9, 6))
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.barplot(data=data[iiCoeff].melt(var_name="Type",
                                        ignore_index=False).reset_index(names="Operation"),
                x='Type', y='value', hue='Operation', width=0.5,
                palette=mycolours[0:3], ax=ax)
    if iiCoeff < len(["Spearman", "Kendall"]) - 1:
        ax.set_xticks(ticks=range(len(labelsLoud)), labels=[])
        ax.set(xlabel=None)
    else:
        ax.set(xlabel="Loudness metric")
        ax.set_xticks(ticks=range(len(labelsLoud)), labels=labelsLoud,
                      rotation=45, ha='right', rotation_mode='anchor')

    ax.set(yticks=np.arange(0, 1.1, 0.1),
           ylabel="Correlation coefficient"
           + "\n" + "(median annoyance rating)")
    ax.set_title("  " + Coeff, y=1.0, pad=-10, loc='left')

    ax.legend(loc='upper right', ncol=len(UASOps))
    ax.grid(alpha=0.15, linestyle='--', axis='y')


# By UAS type
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7.5, 6))
data = [pd.DataFrame(data=[corr_matrixAByUAS[0][1].loc['AnnoyMedian',
                                                       indicesLoudness],
                           corr_matrixAByUAS[1][1].loc['AnnoyMedian',
                                                       indicesLoudness],
                           corr_matrixAByUAS[2][1].loc['AnnoyMedian',
                                                       indicesLoudness]],
                     index=["H520E", "M300", "T150"]),
        pd.DataFrame(data=[corr_matrixAByUAS[0][2].loc['AnnoyMedian',
                                                       indicesLoudness],
                           corr_matrixAByUAS[1][2].loc['AnnoyMedian',
                                                       indicesLoudness],
                           corr_matrixAByUAS[2][2].loc['AnnoyMedian',
                                                       indicesLoudness]],
                     index=["H520E", "M300", "T150"])]
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.heatmap(data=data[iiCoeff], vmin=0.3, vmax=0.9,
                xticklabels=labelsLoud, yticklabels=data[iiCoeff].index,
                ax=ax, cbar_kws=dict(label="Correlation coefficient"
                                     + "\n" + "(median annoyance rating)"))
    ax.set_xticks(ticks=range(len(labelsLoud)), labels=labelsLoud, rotation=30)
    ax.set_title("  " + Coeff, y=1.0, pad=5, loc='left')


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(9, 6))
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.barplot(data=data[iiCoeff].melt(var_name="Type",
                                        ignore_index=False).reset_index(names="UAS Type"),
                x='Type', y='value', hue='UAS Type', width=0.5,
                palette=mycolours[0:3], ax=ax)
    if iiCoeff < len(["Spearman", "Kendall"]) - 1:
        ax.set_xticks(ticks=range(len(labelsLoud)), labels=[])
        ax.set(xlabel=None)
    else:
        ax.set(xlabel="Loudness metric")
        ax.set_xticks(ticks=range(len(labelsLoud)), labels=labelsLoud,
                      rotation=45, ha='right', rotation_mode='anchor')

    ax.set(yticks=np.arange(0, 1.1, 0.1),
           ylabel=("Correlation coefficient"
                   + "\n" + "(median annoyance rating)"))
    ax.set_title("  " + Coeff, y=1.0, pad=-10, loc='left')

    ax.legend(loc='upper right', ncol=len(UASTypes))
    ax.grid(alpha=0.15, linestyle='--', axis='y')

# sound quality

indicesSonQual = ['TonalECMAHMSAvgMaxLR', 'RoughECMAHMS10ExBin',
                  'FluctstrHMS10ExBin', 'SharpAuresISO3AvgMaxLR',
                  'SharpAuresISO3PowAvgMaxLR', 'ImpulsHMSAvgMaxLR',
                  'ImpulsHMSMaxMaxLR', 'ImpulsHMSPowAvgMaxLR']
indicesSonQual = indicesSonQual + (["UAS" + index
                                    for index in indicesSonQual]
                                   + ["Amb" + index
                                      for index in indicesSonQual])

labelsSonQual = ["$T_\\mathrm{mean,ECMA(HMS)}$", "$R_\\mathrm{10,ECMA(HMS)}$",
                 "$F_\\mathrm{10,HMS}$", "$S_\\mathrm{mean,Aures,ISO3}$",
                 "$S_\\mathrm{Aures,ISO3}$", "$I_\\mathrm{mean,HMS}$",
                 "$I_\\mathrm{max,HMS}$", "$I_\\mathrm{HMS}$"]
labelsSonQual = labelsSonQual + (["UAS " + label for label in labelsSonQual]
                                 + ["Amb " + label for label in labelsSonQual])

fig, ax = plt.subplots(figsize=(8, 4))
data = pd.DataFrame(data=[corr_matrixA[1].loc['AnnoyMedian', indicesSonQual],
                          corr_matrixA[2].loc['AnnoyMedian', indicesSonQual]],
                    index=["Spearman", "Kendall"])
sns.barplot(data=data.melt(var_name="Type",
                           ignore_index=False).reset_index(names="Method"),
            x='Type', y='value', hue='Method', width=0.75,
            palette=mycolours[0:2], ax=ax)
ax.set(yticks=np.arange(-0.6, 1.1, 0.1), xlabel="Sound quality metric",
       ylabel="Correlation coefficient" + "\n" + "(median annoyance rating)")
ax.set_xticks(ticks=range(len(labelsSonQual)), labels=labelsSonQual,
              rotation=45, ha='right', rotation_mode='anchor')
ax.legend(loc='best')
ax.grid(alpha=0.15, linestyle='--', axis='y')

# plot UAS sound quality data
labelsSonQualSub = (labelsSonQual[8:11] + [labelsSonQual[12]]
                    + [labelsSonQual[14]])
fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(15, 3.75))
for ii, (index, ax) in enumerate(zip(indicesSonQual[8:11]
                                     + [indicesSonQual[12]]
                                     + [indicesSonQual[14]],
                                     axs.ravel())):
    sns.regplot(data=ptADataForCorr, x=index, y='AnnoyMedian',
                ax=ax, y_jitter=0.1, color=mycolours[ii],
                scatter=True, fit_reg=False, scatter_kws=dict(alpha=0.4))
    if index == 'UASFluctstrHMS10ExBin':
        ax.set(xticks=np.arange(np.round(np.min(ptADataForCorr[index]), 1),
                                np.ceil(np.max(ptADataForCorr[index])*100)/100,
                                np.ceil(np.max(ptADataForCorr[index])*100)/100/5))
    else:
        ax.set(xticks=np.arange(np.round(np.min(ptADataForCorr[index]), 2),
                                np.round(np.ceil(np.max(ptADataForCorr[index])*10)/10, 1),
                                np.round(np.ceil(np.max(ptADataForCorr[index])*10)/10, 1)/5))
    ax.set(yticks=range(0, 11),
           xlabel=labelsSonQualSub[ii],
           ylabel="Median annoyance rating")
    ax.grid(alpha=0.15, linestyle='--')
plt.show()

# plot whole scene sound quality data
labelsSonQualSub = (labelsSonQual[0:3] + [labelsSonQual[4]]
                    + [labelsSonQual[6]])
fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(15, 3.75))
ptADataForCorrSub1 = ptADataForCorr.loc[~ptADataForCorr.index.isin(["A1_CALBIN_Pa.wav",
                                                                    "A2_CALBIN_Pa.wav"])]
ptADataForCorrSub2 = ptADataForCorr.loc[ptADataForCorr.index.isin(["A1_CALBIN_Pa.wav",
                                                                   "A2_CALBIN_Pa.wav"])]
for ii, (index, ax) in enumerate(zip(indicesSonQual[0:3]
                                     + [indicesSonQual[4]]
                                     + [indicesSonQual[6]],
                                     axs.ravel())):
    sns.regplot(data=ptADataForCorrSub1, x=index, y='AnnoyMedian',
                ax=ax, y_jitter=0.1, color=mycolours[ii],
                scatter=True, fit_reg=False, scatter_kws=dict(alpha=0.4))
    sns.regplot(data=ptADataForCorrSub2, x=index, y='AnnoyMedian',
                ax=ax, y_jitter=0.1, color=mycolours[ii], marker='x',
                scatter=True, fit_reg=False, scatter_kws=dict(alpha=0.4))
    if index == 'UASFluctstrHMS10ExBin':
        ax.set(xticks=np.arange(np.round(np.min(ptADataForCorr[index]), 1),
                                np.ceil(np.max(ptADataForCorr[index])*100)/100,
                                np.ceil(np.max(ptADataForCorr[index])*100)/100/5))
    else:
        ax.set(xticks=np.arange(np.round(np.min(ptADataForCorr[index]), 2),
                                np.round(np.ceil(np.max(ptADataForCorr[index])*10)/10, 1),
                                np.round(np.ceil(np.max(ptADataForCorr[index])*10)/10, 1)/5))
    ax.set(yticks=range(0, 11),
           xlabel=labelsSonQualSub[ii],
           ylabel="Median annoyance rating")
    ax.grid(alpha=0.15, linestyle='--')
plt.show()


# By operation
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(13, 6))
data = [pd.DataFrame(data=[corr_matrixAByOp[0][1].loc['AnnoyMedian',
                                                      indicesSonQual],
                           corr_matrixAByOp[1][1].loc['AnnoyMedian',
                                                      indicesSonQual],
                           corr_matrixAByOp[2][1].loc['AnnoyMedian',
                                                      indicesSonQual]],
                     index=["Flyover", "Landing", "Takeoff"]),
        pd.DataFrame(data=[corr_matrixAByOp[0][2].loc['AnnoyMedian',
                                                      indicesSonQual],
                           corr_matrixAByOp[1][2].loc['AnnoyMedian',
                                                      indicesSonQual],
                           corr_matrixAByOp[2][2].loc['AnnoyMedian',
                                                      indicesSonQual]],
                     index=["Flyover", "Landing", "Takeoff"])]
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.heatmap(data=data[iiCoeff], vmin=-0.5, vmax=0.7,
                xticklabels=labelsSonQual, yticklabels=data[iiCoeff].index,
                ax=ax, cmap='seismic',
                cbar_kws=dict(label="Correlation coefficient"
                              + "\n" + "(median annoyance rating)"))
    ax.set_xticks(ticks=range(len(labelsSonQual)), labels=labelsSonQual,
                  rotation=30)
    ax.set_title("  " + Coeff, y=1.0, pad=5, loc='left')


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(13, 6))
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.barplot(data=data[iiCoeff].melt(var_name="Type",
                                        ignore_index=False).reset_index(names="Operation"),
                x='Type', y='value', hue='Operation', width=0.5,
                palette=mycolours[0:3], ax=ax)
    if iiCoeff < len(["Spearman", "Kendall"]) - 1:
        ax.set_xticks(ticks=range(len(labelsSonQual)), labels=[])
        ax.set(xlabel=None)
    else:
        ax.set(xlabel="Sound quality metric")
        ax.set_xticks(ticks=range(len(labelsSonQual)), labels=labelsSonQual,
                      rotation=45, ha='right', rotation_mode='anchor')

    ax.set(yticks=np.arange(-0.6, 1.1, 0.1),
           ylabel="Correlation coefficient"
           + "\n" + "(median annoyance rating)")
    ax.set_title("  " + Coeff, y=1.0, pad=-10, loc='left')

    ax.legend(loc='upper right', ncol=len(UASOps))
    ax.grid(alpha=0.15, linestyle='--', axis='y')


# By UAS type
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(13, 6))
data = [pd.DataFrame(data=[corr_matrixAByUAS[0][1].loc['AnnoyMedian',
                                                       indicesSonQual],
                           corr_matrixAByUAS[1][1].loc['AnnoyMedian',
                                                       indicesSonQual],
                           corr_matrixAByUAS[2][1].loc['AnnoyMedian',
                                                       indicesSonQual]],
                     index=["H520E", "M300", "T150"]),
        pd.DataFrame(data=[corr_matrixAByUAS[0][2].loc['AnnoyMedian',
                                                       indicesSonQual],
                           corr_matrixAByUAS[1][2].loc['AnnoyMedian',
                                                       indicesSonQual],
                           corr_matrixAByUAS[2][2].loc['AnnoyMedian',
                                                       indicesSonQual]],
                     index=["H520E", "M300", "T150"])]
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.heatmap(data=data[iiCoeff], vmin=-0.5, vmax=0.8,
                xticklabels=labelsSonQual, yticklabels=data[iiCoeff].index,
                ax=ax, cmap='seismic',
                cbar_kws=dict(label="Correlation coefficient"
                              + "\n" + "(median annoyance rating)"))
    ax.set_xticks(ticks=range(len(labelsSonQual)), labels=labelsSonQual,
                  rotation=30)
    ax.set_title("  " + Coeff, y=1.0, pad=5, loc='left')


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(13, 6))
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.barplot(data=data[iiCoeff].melt(var_name="Type",
                                        ignore_index=False).reset_index(names="UAS Type"),
                x='Type', y='value', hue='UAS Type', width=0.5,
                palette=mycolours[0:3], ax=ax)
    if iiCoeff < len(["Spearman", "Kendall"]) - 1:
        ax.set_xticks(ticks=range(len(labelsSonQual)), labels=[])
        ax.set(xlabel=None)
    else:
        ax.set(xlabel="Sound quality metric")
        ax.set_xticks(ticks=range(len(labelsSonQual)), labels=labelsSonQual,
                      rotation=45, ha='right', rotation_mode='anchor')

    ax.set(yticks=np.arange(-0.6, 1.1, 0.1),
           ylabel="Correlation coefficient"
           + "\n" + "(median annoyance rating)")
    ax.set_title("  " + Coeff, y=1.0, pad=-10, loc='left')

    ax.legend(loc='upper right', ncol=len(UASTypes))
    ax.grid(alpha=0.15, linestyle='--', axis='y')


# -----------------------------------------------------------------------------
# Part B
# -----------------------------------------------------------------------------

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

# Add jittered dots
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

# Add jittered dots
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
# Add jittered dots
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

# scatter plots and correlation analysis
# --------------------------------------
# reduce dataset
includecols = ["LAeqMaxLR", "LAEMaxLR", "LAFmaxMaxLR",
               "LAF5ExMaxLR", "LAF10ExMaxLR", "LAF25ExMaxLR", "LAF50ExMaxLR",
               "LAF75ExMaxLR", "LAF90ExMaxLR", "LAF95ExMaxLR", "LASmaxMaxLR",
               "EPNLMaxLR", "EPNLNoTMaxLR", "PNLTmaxMaxLR",
               "LoudECMAHMSPowAvgBin", "TonalECMAHMSAvgMaxLR",
               "RoughECMAHMS10ExBin", "FluctstrHMS10ExBin",
               "LoudISO105ExMaxLR", "LoudISO1PowAvgMaxLR",
               "LoudISO3PowAvgMaxLR", "SharpAuresISO3AvgMaxLR",
               "SharpAuresISO3PowAvgMaxLR", "ImpulsHMSPowAvgMaxLR",
               "ImpulsHMSAvgMaxLR", "ImpulsHMSMaxMaxLR"]
includecols = (includecols +
               ["UAS" + includecol for includecol in includecols])
includecols = includecols + ["SNRlevel", "UASLAeqdiffAmbLAF90",
                             "UASLASmaxdiffAmbLAF90", "UASLASmaxdiffAmbLAF50",
                             "UASLASmaxdiffAmbLAeq", "UASPartLoudGMSTPowAvg",
                             "EPNLdiff", "EPNLNoTdiff", "PNLTmaxdiff"]
includecols = includecols + ['UASEvents']
includecols = (["AnnoyMedian", "AnnoyMean", "ValenceMedian", "ValenceMean",
                "ArousalMedian", "ArousalMean"]
               + includecols)

ptBDataForCorr = dataByStimTestB.loc[:, includecols]
# convert categorical data to numerical
ptBDataForCorr['SNRlevel'] = ptBDataForCorr['SNRlevel'].cat.codes
ptBDataForCorr.loc[ptBDataForCorr.loc[:, 'SNRlevel']
                   == -1, 'SNRlevel'] = np.nan
ptBDataForCorr['SNRlevel'] = ptBDataForCorr['SNRlevel'].astype(float)

UASTypes = dataByStimTestB['UASType'].unique()[1:]
ptBDataForCorrByUAS = [dataTemp for dataTemp in
                       [dataByStimTestB.loc[dataByStimTestB['UASType'] == UAS,
                                            includecols] for UAS
                        in UASTypes]]



# run correlation analysis
corr_matrixB = [ptBDataForCorr.corr(method=corr_method, numeric_only=False)
                for corr_method in ['pearson', 'spearman', 'kendall']]
corr_matrixBByUAS = [[dataByUAS.corr(method=corr_method,
                                     numeric_only=False)
                      for corr_method in ['pearson', 'spearman', 'kendall']]
                     for dataByUAS in ptBDataForCorrByUAS]

# re-run correlation analysis to get p-values
corr_ListB_pr = pd.DataFrame()
corr_ListB_sp = pd.DataFrame()
corr_ListB_kd = pd.DataFrame()
feat1s = [None]*len(ptBDataForCorr.columns)**2
feat2s = [None]*len(ptBDataForCorr.columns)**2
corrs_pr = np.zeros(len(ptBDataForCorr.columns)**2)
corrs_sp = np.zeros(len(ptBDataForCorr.columns)**2)
corrs_kd = np.zeros(len(ptBDataForCorr.columns)**2)
pvalues_pr = np.zeros(len(ptBDataForCorr.columns)**2)
pvalues_sp = np.zeros(len(ptBDataForCorr.columns)**2)
pvalues_kd = np.zeros(len(ptBDataForCorr.columns)**2)

count = 0
for feat1 in ptBDataForCorr.columns:
    for feat2 in ptBDataForCorr.columns:
        if feat1 == feat2:
            corr_pr, pvalue_pr = 1, 0
            corr_kd, pvalue_kd = 1, 0
        else:
            # avoid nan issues and speed up Kendall's tau
            mask = ~np.logical_or(ptBDataForCorr[feat1].isna(),
                                  ptBDataForCorr[feat2].isna())

            corr_pr, pvalue_pr = stats.pearsonr(ptBDataForCorr.loc[mask,
                                                                   feat1].astype(float),
                                                ptBDataForCorr.loc[mask,
                                                                   feat2].astype(float))
            corr_kd, pvalue_kd = stats.kendalltau(ptBDataForCorr.loc[mask,
                                                                     feat1].astype(float),
                                                  ptBDataForCorr.loc[mask,
                                                                     feat2].astype(float))
        corr_sp, pvalue_sp = stats.spearmanr(ptBDataForCorr[feat1].astype(float),
                                             ptBDataForCorr[feat2].astype(float),
                                             nan_policy='omit')

        feat1s[count] = feat1
        feat2s[count] = feat2
        corrs_pr[count] = corr_pr
        pvalues_pr[count] = pvalue_pr
        corrs_sp[count] = corr_sp
        pvalues_sp[count] = pvalue_sp
        corrs_kd[count] = corr_kd
        pvalues_kd[count] = pvalue_kd

        count += 1

# end of for loop for correlation analysis

corr_ListB_pr['Variable_1'] = feat1s
corr_ListB_pr['Variable_2'] = feat2s
corr_ListB_pr['Correlation'] = corrs_pr
corr_ListB_pr['p_value'] = pvalues_pr

corr_ListB_sp['Variable_1'] = feat1s
corr_ListB_sp['Variable_2'] = feat2s
corr_ListB_sp['Correlation'] = corrs_sp
corr_ListB_sp['p_value'] = pvalues_sp

corr_ListB_kd['Variable_1'] = feat1s
corr_ListB_kd['Variable_2'] = feat2s
corr_ListB_kd['Correlation'] = corrs_kd
corr_ListB_kd['p_value'] = pvalues_kd

corr_matrixB_p = [corrlist.pivot(columns='Variable_2',
                                 index='Variable_1',
                                 values='p_value')[corr_matrixB[1].columns].reindex(corr_matrixB[1].index)
                  for corrlist in [corr_ListB_pr,
                                   corr_ListB_sp,
                                   corr_ListB_kd]]

# response correlation analysis
# -----------------------------

# plot annoyance correlations
# sound level indices
indicesAcoustic = ['LAeqMaxLR', 'LAFmaxMaxLR', 'LAF5ExMaxLR', 'LAF10ExMaxLR',
                   'LAF25ExMaxLR', 'LAF50ExMaxLR', 'LAF75ExMaxLR',
                   'LAF90ExMaxLR', 'LAF95ExMaxLR', 'LASmaxMaxLR', 'EPNLMaxLR',
                   'EPNLNoTMaxLR', 'PNLTmaxMaxLR']
indicesAcoustic = indicesAcoustic + (["UAS" + index
                                      for index in indicesAcoustic])

labelsAcoustic = ["$L_\\mathrm{Aeq}$", "$L_\\mathrm{AFmax}$",
                  "$L_\\mathrm{AF5}$", "$L_\\mathrm{AF10}$",
                  "$L_\\mathrm{AF25}$", "$L_\\mathrm{AF50}$",
                  "$L_\\mathrm{AF75}$", "$L_\\mathrm{AF90}$",
                  "$L_\\mathrm{AF95}$", "$L_\\mathrm{ASmax}$",
                  "EPNL", "EPNL,$_\\mathrm{no\\ tone}$",
                  "PNLT$_\\mathrm{max}$"]
labelsAcoustic = labelsAcoustic + (["UAS " + label
                                    for label in labelsAcoustic])

data = pd.DataFrame(data=[corr_matrixB[1].loc['AnnoyMedian', indicesAcoustic],
                          corr_matrixB[2].loc['AnnoyMedian', indicesAcoustic]],
                    index=["Spearman", "Kendall"])

fig, ax = plt.subplots(figsize=(7, 4))
sns.barplot(data=data.melt(var_name="Type",
                           ignore_index=False).reset_index(names="Method"),
            x='Type', y='value', hue='Method', width=0.75,
            palette=mycolours[0:2], ax=ax)
ax.set(yticks=np.arange(-0.1, 1.1, 0.1), xlabel="Sound level index",
       ylabel="Correlation coefficient" + "\n" + "(median annoyance)")
ax.set_xticks(ticks=range(len(labelsAcoustic)), labels=labelsAcoustic,
              rotation=45, ha='right', rotation_mode='anchor')
ax.legend(loc='best')
ax.grid(alpha=0.15, linestyle='--', axis='y')
plt.show()

fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(14, 7))
labelsAcousticSub = ([labelsAcoustic[13]]
                     + labelsAcoustic[16:22] + [labelsAcoustic[23]])
for ii, (index, ax) in enumerate(zip([indicesAcoustic[13]]
                                     + indicesAcoustic[16:22]
                                     + [indicesAcoustic[23]],
                                     axs.ravel())):
    sns.regplot(data=ptBDataForCorr, x=index, y='AnnoyMedian',
                ax=ax, y_jitter=0.1, color=mycolours[ii],
                scatter=True, fit_reg=False, scatter_kws=dict(alpha=0.4))
    ax.set(xticks=range(int(np.min(ptBDataForCorr[index])/5)*5,
                        round(np.ceil(np.max(ptBDataForCorr[index]))/5)*5 + 5,
                        5), yticks=range(0, 11),
           xlabel=labelsAcousticSub[ii] + ", dB",
           ylabel="Median annoyance rating")
    ax.grid(alpha=0.15, linestyle='--')


# By UAS type
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(10, 6))
UASTypes = dataByStimTestB['UASType'].unique()[1:]
data = [pd.DataFrame(data=[corr_matrixBByUAS[0][1].loc['AnnoyMedian',
                                                       indicesAcoustic],
                           corr_matrixBByUAS[1][1].loc['AnnoyMedian',
                                                       indicesAcoustic]],
                     index=["H520E", "T150"]),
        pd.DataFrame(data=[corr_matrixBByUAS[0][2].loc['AnnoyMedian',
                                                       indicesAcoustic],
                           corr_matrixBByUAS[1][2].loc['AnnoyMedian',
                                                       indicesAcoustic]],
                     index=["H520E", "T150"])]
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.heatmap(data=data[iiCoeff], vmin=-0.4, vmax=0.9,
                xticklabels=labelsAcoustic, yticklabels=data[iiCoeff].index,
                ax=ax, cmap='seismic',
                cbar_kws=dict(label="Correlation coefficient"
                              + "\n" + "(median annoyance rating)"))
    ax.set_xticks(ticks=range(len(labelsAcoustic)), labels=labelsAcoustic,
                  rotation=30, ha='right', rotation_mode='anchor')
    ax.set_title("  " + Coeff, y=1.0, pad=5, loc='left')


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(12, 6))
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.barplot(data=data[iiCoeff].melt(var_name="Type",
                                        ignore_index=False).reset_index(names="UAS Type"),
                x='Type', y='value', hue='UAS Type', width=0.5,
                palette=mycolours[0:2], ax=ax)
    if iiCoeff < len(["Spearman", "Kendall"]) - 1:
        ax.set_xticks(ticks=range(len(labelsAcoustic)), labels=[])
        ax.set(xlabel=None)
    else:
        ax.set(xlabel="UAS sound vs ambient sound")
        ax.set_xticks(ticks=range(len(labelsAcoustic)), labels=labelsAcoustic,
                      rotation=45, ha='right', rotation_mode='anchor')

    ax.set(yticks=np.arange(-0.5, 1.1, 0.1),
           ylabel="Correlation coefficient"
           + "\n" + "(median annoyance rating)")
    ax.set_title("  " + Coeff, y=1.0, pad=-14, loc='left')

    ax.legend(loc='upper right', ncol=len(UASTypes))
    ax.grid(alpha=0.15, linestyle='--', axis='y')


# loudness metrics
indicesLoudness = ['LoudECMAHMSPowAvgBin', 'LoudISO105ExMaxLR',
                   'LoudISO1PowAvgMaxLR', 'LoudISO3PowAvgMaxLR']
indicesLoudness = indicesLoudness + (["UAS" + index
                                      for index in indicesLoudness])

labelsLoud = ["$N_\\mathrm{ECMA(HMS)}$", "$N_{5\\mathrm{,ISO1}}$",
              "$N_\\mathrm{ISO1}$", "$N_\\mathrm{ISO3}$"]
labelsLoud = labelsLoud + (["UAS " + label for label in labelsLoud])

fig, ax = plt.subplots(figsize=(5, 4))
data = pd.DataFrame(data=[corr_matrixB[1].loc['AnnoyMedian', indicesLoudness],
                          corr_matrixB[2].loc['AnnoyMedian', indicesLoudness]],
                    index=["Spearman", "Kendall"])
sns.barplot(data=data.melt(var_name="Type",
                           ignore_index=False).reset_index(names="Method"),
            x='Type', y='value', hue='Method', width=0.75,
            palette=mycolours[0:2], ax=ax)
ax.set(yticks=np.arange(0, 1.1, 0.1), xlabel="Loudness metric",
       ylabel="Correlation coefficient" + "\n" + "(median annoyance)")
ax.set_xticks(ticks=range(len(labelsLoud)), labels=labelsLoud, rotation=45,
              ha='right', rotation_mode='anchor')
ax.legend(loc='best', ncol=2)
ax.grid(alpha=0.15, linestyle='--', axis='y')


labelsLoudSub = [labelsLoud[4]] + [labelsLoud[7]]
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(7, 3.75))
for ii, (index, ax) in enumerate(zip([indicesLoudness[4]]
                                     + [indicesLoudness[7]],
                                     axs.ravel())):
    sns.regplot(data=ptBDataForCorr, x=index, y='AnnoyMedian',
                ax=ax, y_jitter=0.1, color=mycolours[ii],
                scatter=True, fit_reg=False, scatter_kws=dict(alpha=0.4))
    ax.set(xticks=range(int(np.min(ptBDataForCorr[index])/2)*2,
                        round(np.ceil(np.max(ptBDataForCorr[index]))/2)*2 + 2,
                        2), yticks=range(0, 11),
           xlabel=labelsLoudSub[ii] + ", sones",
           ylabel="Median annoyance rating")
    ax.grid(alpha=0.15, linestyle='--')


# By UAS type
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(7.5, 6))
UASTypes = dataByStimTestB['UASType'].unique()[1:]
data = [pd.DataFrame(data=[corr_matrixBByUAS[0][1].loc['AnnoyMedian',
                                                       indicesLoudness],
                           corr_matrixBByUAS[1][1].loc['AnnoyMedian',
                                                       indicesLoudness]],
                     index=["H520E", "T150"]),
        pd.DataFrame(data=[corr_matrixBByUAS[0][2].loc['AnnoyMedian',
                                                       indicesLoudness],
                           corr_matrixBByUAS[1][2].loc['AnnoyMedian',
                                                       indicesLoudness]],
                     index=["H520E", "T150"])]
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.heatmap(data=data[iiCoeff], vmin=0.3, vmax=0.9,
                xticklabels=labelsLoud, yticklabels=data[iiCoeff].index,
                ax=ax, cbar_kws=dict(label="Correlation coefficient"
                                     + "\n" + "(median annoyance rating)"))
    ax.set_xticks(ticks=range(len(labelsLoud)), labels=labelsLoud, rotation=30)
    ax.set_title("  " + Coeff, y=1.0, pad=5, loc='left')


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.barplot(data=data[iiCoeff].melt(var_name="Type",
                                        ignore_index=False).reset_index(names="UAS Type"),
                x='Type', y='value', hue='UAS Type', width=0.5,
                palette=mycolours[0:2], ax=ax)
    if iiCoeff < len(["Spearman", "Kendall"]) - 1:
        ax.set_xticks(ticks=range(len(labelsLoud)), labels=[])
        ax.set(xlabel=None)
    else:
        ax.set(xlabel="UAS sound vs ambient sound")
        ax.set_xticks(ticks=range(len(labelsLoud)), labels=labelsLoud,
                      rotation=45, ha='right', rotation_mode='anchor')

    ax.set(yticks=np.arange(0, 1.1, 0.1),
           ylabel="Correlation coefficient" + "\n" + "(median annoyance rating)")
    ax.set_title("  " + Coeff, y=1.0, pad=-10, loc='left')

    ax.legend(loc='upper right', ncol=len(UASTypes))
    ax.grid(alpha=0.15, linestyle='--', axis='y')


# sound quality

indicesSonQual = ['TonalECMAHMSAvgMaxLR', 'RoughECMAHMS10ExBin',
                  'FluctstrHMS10ExBin', 'SharpAuresISO3AvgMaxLR',
                  'SharpAuresISO3PowAvgMaxLR', 'ImpulsHMSAvgMaxLR',
                  'ImpulsHMSMaxMaxLR', 'ImpulsHMSPowAvgMaxLR']
indicesSonQual = indicesSonQual + (["UAS" + index
                                    for index in indicesSonQual])

labelsSonQual = ["$T_\\mathrm{mean,ECMA(HMS)}$", "$R_\\mathrm{10,ECMA(HMS)}$",
                 "$F_\\mathrm{10,HMS}$", "$S_\\mathrm{mean,Aures,ISO3}$",
                 "$S_\\mathrm{Aures,ISO3}$", "$I_\\mathrm{mean,HMS}$",
                 "$I_\\mathrm{max,HMS}$", "$I_\\mathrm{HMS}$"]
labelsSonQual = labelsSonQual + (["UAS " + label for label in labelsSonQual])

fig, ax = plt.subplots(figsize=(6, 4))
data = pd.DataFrame(data=[corr_matrixB[1].loc['AnnoyMedian', indicesSonQual],
                          corr_matrixB[2].loc['AnnoyMedian', indicesSonQual]],
                    index=["Spearman", "Kendall"])
sns.barplot(data=data.melt(var_name="Type",
                           ignore_index=False).reset_index(names="Method"),
            x='Type', y='value', hue='Method', width=0.75,
            palette=mycolours[0:2], ax=ax)
ax.set(yticks=np.arange(-0.8, 1.1, 0.1), xlabel="Sound quality metric",
       ylabel="Correlation coefficient" + "\n" + "(median annoyance rating)")
ax.set_xticks(ticks=range(len(labelsSonQual)), labels=labelsSonQual,
              rotation=45, ha='right', rotation_mode='anchor')
ax.legend(loc='best')
ax.grid(alpha=0.15, linestyle='--', axis='y')
plt.show()

# plot UAS sound quality data
labelsSonQualSub = (labelsSonQual[8:11] + [labelsSonQual[12]]
                    + [labelsSonQual[14]])
fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(15, 3.75))
for ii, (index, ax) in enumerate(zip(indicesSonQual[8:11]
                                     + [indicesSonQual[12]]
                                     + [indicesSonQual[14]],
                                     axs.ravel())):
    sns.regplot(data=ptBDataForCorr, x=index, y='AnnoyMedian',
                ax=ax, y_jitter=0.1, color=mycolours[ii],
                scatter=True, fit_reg=False, scatter_kws=dict(alpha=0.4))
    if index == 'UASFluctstrHMS10ExBin':
        ax.set(xticks=np.arange(np.round(np.min(ptBDataForCorr[index]), 1),
                                np.ceil(np.max(ptBDataForCorr[index])*100)/100,
                                np.ceil(np.max(ptBDataForCorr[index])*100)/100/5))
    else:
        ax.set(xticks=np.arange(np.round(np.min(ptBDataForCorr[index]), 2),
                                np.round(np.ceil(np.max(ptBDataForCorr[index])*10)/10, 1),
                                np.round(np.ceil(np.max(ptBDataForCorr[index])*10)/10, 1)/5))
    ax.set(yticks=range(0, 11),
           xlabel=labelsSonQualSub[ii],
           ylabel="Median annoyance rating")
    ax.grid(alpha=0.15, linestyle='--')
plt.show()

# plot whole scene sound quality data
labelsSonQualSub = (labelsSonQual[0:3] + [labelsSonQual[4]]
                    + [labelsSonQual[6]])
fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(15, 3.75))
for ii, (index, ax) in enumerate(zip(indicesSonQual[0:3]
                                     + [indicesSonQual[4]]
                                     + [indicesSonQual[6]],
                                     axs.ravel())):
    sns.regplot(data=ptBDataForCorr, x=index, y='AnnoyMedian',
                ax=ax, y_jitter=0.1, color=mycolours[ii],
                scatter=True, fit_reg=False, scatter_kws=dict(alpha=0.4))
    if index == 'UASFluctstrHMS10ExBin':
        ax.set(xticks=np.arange(np.round(np.min(ptBDataForCorr[index]), 1),
                                np.ceil(np.max(ptBDataForCorr[index])*100)/100,
                                np.ceil(np.max(ptBDataForCorr[index])*100)/100/5))
    else:
        ax.set(xticks=np.arange(np.round(np.min(ptBDataForCorr[index]), 2),
                                np.round(np.ceil(np.max(ptBDataForCorr[index])*10)/10, 1),
                                np.round(np.ceil(np.max(ptBDataForCorr[index])*10)/10, 1)/5))
    ax.set(yticks=range(0, 11),
           xlabel=labelsSonQualSub[ii],
           ylabel="Median annoyance rating")
    ax.grid(alpha=0.15, linestyle='--')
plt.show()

# By UAS type
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(13, 6))
UASTypes = dataByStimTestB['UASType'].unique()[1:]
data = [pd.DataFrame(data=[corr_matrixBByUAS[0][1].loc['AnnoyMedian',
                                                       indicesSonQual],
                           corr_matrixBByUAS[1][1].loc['AnnoyMedian',
                                                       indicesSonQual],
                           corr_matrixBByUAS[2][1].loc['AnnoyMedian',
                                                       indicesSonQual]],
                     index=["H520E", "M300", "T150"]),
        pd.DataFrame(data=[corr_matrixBByUAS[0][2].loc['AnnoyMedian',
                                                       indicesSonQual],
                           corr_matrixBByUAS[1][2].loc['AnnoyMedian',
                                                       indicesSonQual],
                           corr_matrixBByUAS[2][2].loc['AnnoyMedian',
                                                       indicesSonQual]],
                     index=["H520E", "M300", "T150"])]
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.heatmap(data=data[iiCoeff], vmin=-0.5, vmax=0.8,
                xticklabels=labelsSonQual, yticklabels=data[iiCoeff].index,
                ax=ax, cmap='seismic',
                cbar_kws=dict(label="Correlation coefficient"
                              + "\n" + "(median annoyance rating)"))
    ax.set_xticks(ticks=range(len(labelsSonQual)), labels=labelsSonQual,
                  rotation=30)
    ax.set_title("  " + Coeff, y=1.0, pad=5, loc='left')


fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(13, 6))
for iiCoeff, (Coeff, ax) in enumerate(zip(["Spearman", "Kendall"],
                                          axs.ravel())):
    sns.barplot(data=data[iiCoeff].melt(var_name="Type",
                                        ignore_index=False).reset_index(names="UAS Type"),
                x='Type', y='value', hue='UAS Type', width=0.5,
                palette=mycolours[0:3], ax=ax)
    if iiCoeff < len(["Spearman", "Kendall"]) - 1:
        ax.set_xticks(ticks=range(len(labelsSonQual)), labels=[])
        ax.set(xlabel=None)
    else:
        ax.set(xlabel="Sound quality metric")
        ax.set_xticks(ticks=range(len(labelsSonQual)), labels=labelsSonQual,
                      rotation=45, ha='right', rotation_mode='anchor')

    ax.set(yticks=np.arange(-0.6, 1.1, 0.1),
           ylabel="Correlation coefficient"
           + "\n" + "(median annoyance rating)")
    ax.set_title("  " + Coeff, y=1.0, pad=-10, loc='left')

    ax.legend(loc='upper right', ncol=len(UASTypes))
    ax.grid(alpha=0.15, linestyle='--', axis='y')
