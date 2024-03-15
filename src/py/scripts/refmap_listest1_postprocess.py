# -*- coding: utf-8 -*-

# script


# -----
# Setup
# -----

# import statements
import sys
import os
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import QFileDialog, QApplication
import librosa
import dsp.filterFuncs
from scipy import stats

# enable copy-on-write mode for Pandas (will be default from Pandas 3.0)
pd.options.mode.copy_on_write = True

# output variables
indicesAcoustic = ["LAeqMaxLR", "LAEMaxLR", "LAFmaxMaxLR",
                   "LAF5ExMaxLR", "LAF10ExMaxLR", "LAF25ExMaxLR",
                   "LAF50ExMaxLR", "LAF75ExMaxLR", "LAF90ExMaxLR",
                   "LAF95ExMaxLR", "LASmaxMaxLR"]

# skip signal start and end for time-aggregation
start_skipT = 0.5
end_skipT = 0.5

# ------------------------------------------
# Acoustic metrics single values calculation
# ------------------------------------------

# open wav file selection dialog and assign filepaths to list
# PROJECT NOTE: the calibrated files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Calibration/Post_calib_recs/Take2/Pa_calib?csf=1&web=1&e=oGssZv
app = QApplication(sys.argv)
fileExts = "*.wav"
filelist = list(QFileDialog.getOpenFileNames(filter=fileExts))[0]
filelist.sort()
filenames = [filepath.split('/')[-1] for filepath in filelist]

dataByStim = pd.DataFrame(index=filenames, columns=indicesAcoustic, dtype=float)

# loop over files to analyse
for ii, file in enumerate(filelist):
    # load recording
    signal, sampleRatein = librosa.load(file, sr=None, mono=False)
    signal = np.transpose(signal)

    # apply weighting filters
    signalA = dsp.filterFuncs.A_weight_T(signal, sampleRatein)
    signalmagAF = dsp.filterFuncs.time_weight(signalA, sampleRatein, tau=0.125)
    signalmagAS = dsp.filterFuncs.time_weight(signalA, sampleRatein, tau=1)

    # skip samples
    start_skips = int(start_skipT*sampleRatein)
    end_skips = int(end_skipT*sampleRatein)

    # calculate energy metrics
    signalLAeq = 20*np.log10(np.sqrt((signalA[start_skips:
                                              -end_skips]**2).mean(axis=0))/2e-5)
    signalLAE = (signalLAeq
                 + 10*np.log10(len(signalA[start_skips:
                                           -end_skips])/sampleRatein))

    # calculate percentile metrics
    signaldBAF = 20*np.log10(signalmagAF[start_skips:-end_skips]/2e-5)
    signaldBAS = 20*np.log10(signalmagAS[start_skips:-end_skips]/2e-5)
    signalLAFmax = signaldBAF.max(axis=0)
    signalLAF5 = np.percentile(signaldBAF, q=95, axis=0)
    signalLAF10 = np.percentile(signaldBAF, q=90, axis=0)
    signalLAF25 = np.percentile(signaldBAF, q=75, axis=0)
    signalLAF50 = np.percentile(signaldBAF, q=50, axis=0)
    signalLAF75 = np.percentile(signaldBAF, q=25, axis=0)
    signalLAF90 = np.percentile(signaldBAF, q=10, axis=0)
    signalLAF95 = np.percentile(signaldBAF, q=5, axis=0)
    signalLASmax = signaldBAS.max(axis=0)

    # take max L/R ear only
    dataByStim.loc[filenames[ii], 'LAeqMaxLR'] = signalLAeq.max()
    dataByStim.loc[filenames[ii], 'LAEMaxLR'] = signalLAE.max()
    dataByStim.loc[filenames[ii], 'LAFmaxMaxLR'] = signalLAFmax.max()
    dataByStim.loc[filenames[ii], 'LAF5ExMaxLR'] = signalLAF5.max()
    dataByStim.loc[filenames[ii], 'LAF10ExMaxLR'] = signalLAF10.max()
    dataByStim.loc[filenames[ii], 'LAF25ExMaxLR'] = signalLAF25.max()
    dataByStim.loc[filenames[ii], 'LAF50ExMaxLR'] = signalLAF50.max()
    dataByStim.loc[filenames[ii], 'LAF75ExMaxLR'] = signalLAF75.max()
    dataByStim.loc[filenames[ii], 'LAF90ExMaxLR'] = signalLAF90.max()
    dataByStim.loc[filenames[ii], 'LAF95ExMaxLR'] = signalLAF95.max()
    dataByStim.loc[filenames[ii], 'LASmaxMaxLR'] = signalLASmax.max()

# end of for loop over signal wav files


# ------------------------------------------------
# Psychoacoustic metrics single values calculation
# ------------------------------------------------

# output variables
indicesPsycho = ["LoudECMAHMSPowAvgBin",
                 "TonalECMAHMSAvgMaxLR",
                 "RoughECMAHMS10ExBin",
                 "FluctstrHMS10ExBin",
                 "LoudISO105ExMaxLR",
                 "LoudISO1PowAvgMaxLR",
                 "LoudISO3PowAvgMaxLR",
                 "SharpAuresISO3AvgMaxLR",
                 "SharpAuresISO3PowAvgMaxLR",
                 "ImpulsHMSPowAvgMaxLR",
                 "ImpulsHMSAvgMaxLR",
                 "ImpulsHMSMaxMaxLR"]

dataByStim = pd.concat([dataByStim, pd.DataFrame(index=dataByStim.index,
                                                 columns=indicesPsycho,
                                                 dtype=float)], axis=1)

# open xlsx file selection dialog and assign filepaths to list
# PROJECT NOTE: the calculated files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Analysis/ArtemiS/Output?csf=1&web=1&e=0TAPVz
app = QApplication(sys.argv)
fileExts = "*.xlsx"
filelist = list(QFileDialog.getOpenFileNames(filter=fileExts))[0]
filelist.sort()

# loop over files to analyse
for file in filelist:
    print(file.split('/')[-1])
    workbookdata = pd.read_excel(io=file, sheet_name=None)
    stimulus = workbookdata['Sheet1'].columns[5].split(sep='\'')[1]
    # Calculate ECMA-418-2:2022 Sottek Hearing Model overall loudness from
    # 2-channel specific loudness
    # left channel
    specificLoudnessHMSL = pd.DataFrame(workbookdata['Sheet1'].iloc[14:, 1:].values,
                                        columns=workbookdata['Sheet1'].iloc[13, 1:],
                                        index=workbookdata['Sheet1'].iloc[14:, 0])
    # right channel
    specificLoudnessHMSR = pd.DataFrame(workbookdata['Sheet2'].iloc[14:, 1:].values,
                                        columns=workbookdata['Sheet2'].iloc[13, 1:],
                                        index=workbookdata['Sheet2'].iloc[14:, 0])
    # binaural specific loudness (ECMA-418-2:2022 Equation 118)
    specificLoudnessHMSBin = ((specificLoudnessHMSL**2
                               + specificLoudnessHMSR**2)/2).pow(1/2)
    # binaural time-dependent loudness (ECMA-418-2:2022 Equation 116)
    loudnessHMSTimeVar = specificLoudnessHMSBin.sum(axis=0)*0.5
    # binaural overall (power-averaged) loudness (ECMA-418-2:2022 Equation 117)
    loudnessHMSPowAvg = ((loudnessHMSTimeVar.iloc[int(np.ceil(187.5*start_skipT)):
                                                  -int(np.ceil(187.5*end_skipT))]**(1/np.log10(2))).sum()
                         / len(loudnessHMSTimeVar.iloc[int(np.ceil(187.5*start_skipT)):
                                                       -int(np.ceil(187.5*end_skipT))]))**np.log10(2)

    # Calculate ECMA-418-2:2022 Sottek Hearing Model overall tonality from
    # 2-channel specific tonality
    # left channel
    specificTonalityHMSL = pd.DataFrame(workbookdata['Sheet3'].iloc[14:, 1:].values,
                                        columns=workbookdata['Sheet3'].iloc[13, 1:],
                                        index=workbookdata['Sheet3'].iloc[14:, 0])
    # right channel
    specificTonalityHMSR = pd.DataFrame(workbookdata['Sheet4'].iloc[14:, 1:].values,
                                        columns=workbookdata['Sheet4'].iloc[13, 1:],
                                        index=workbookdata['Sheet4'].iloc[14:, 0])
    # 2-channel time-varing tonality (max, not integration)
    tonalityHMSTimeVar = pd.concat([specificTonalityHMSL.max(axis=0),
                                    specificTonalityHMSR.max(axis=0)],
                                   axis=1)
    # 2-channel time-averaged tonality (omitting T<=0.02)
    tonalityHMSAvgL = tonalityHMSTimeVar.iloc[int(np.ceil(187.5*start_skipT)):
                                              -int(np.ceil(187.5*end_skipT)),
                                              0][tonalityHMSTimeVar.iloc[int(np.ceil(187.5*start_skipT)):
                                                                         -int(np.ceil(187.5*end_skipT)), 0]
                                                                         >= 0.02].mean(axis=0)
    tonalityHMSAvgR = tonalityHMSTimeVar.iloc[int(np.ceil(187.5*start_skipT)):
                                              -int(np.ceil(187.5*end_skipT)),
                                              1][tonalityHMSTimeVar.iloc[int(np.ceil(187.5*start_skipT)):
                                                                         -int(np.ceil(187.5*end_skipT)), 1]
                                                                         >= 0.02].mean(axis=0)
    tonalityHMSAvgMaxLR = max(tonalityHMSAvgL, tonalityHMSAvgR)

    # Calculate ECMA-418-2:2022 Sottek Hearing Model binaural overall roughness
    # from binaural specific roughness
    specificRoughnessHMSBin = pd.DataFrame(workbookdata['Sheet5'].iloc[14:, 1:].values,
                                        columns=workbookdata['Sheet5'].iloc[13, 1:],
                                        index=workbookdata['Sheet5'].iloc[14:, 0])
    # binaural time-varying roughness
    roughnessHMSTimeVar = specificRoughnessHMSBin.sum(axis=0)*0.5
    # binaural overall (90th percentile = 10% exceeded) roughness
    roughnessHMS10Ex = roughnessHMSTimeVar.iloc[int(np.ceil(50*start_skipT)):
                                                -int(np.ceil(50*end_skipT))].quantile(q=0.9)

    # Calculate overall Sottek Hearing Model fluctuation strength from
    # 2-channel specific fluctuation strength 
    # calculate according to ECMA-418-2:2022 approach for roughness
    # left channel
    specificFluctStrHMSL = pd.DataFrame(workbookdata['Sheet6'].iloc[14:, 1:].values,
                                        columns=workbookdata['Sheet6'].iloc[13, 1:],
                                        index=workbookdata['Sheet6'].iloc[14:, 0])
    # right channel
    specificFluctStrHMSR = pd.DataFrame(workbookdata['Sheet7'].iloc[14:, 1:].values,
                                        columns=workbookdata['Sheet7'].iloc[13, 1:],
                                        index=workbookdata['Sheet7'].iloc[14:, 0])
    # binaural specific fluctuation strength
    # (using ECMA-418-2:2022 Equation 112 for roughness)
    specificFluctStrHMSBin = ((specificFluctStrHMSL**2
                               + specificFluctStrHMSR**2)/2).pow(1/2)
    # binaural time-dependent fluctuation strength
    # (using ECMA-418-2:2022 Equation 111 for roughness)
    fluctStrHMSTimeVar = specificFluctStrHMSBin.sum(axis=0)*0.5
    # binaural overall (90th percentile = 10% exceeded) fluctuation strength
    # (using ECMA-418-2:2022 Section 7.1.8 for roughness)
    fluctStrHMSTime10Ex = fluctStrHMSTimeVar.iloc[int(np.ceil(start_skipT*(229390681/4e6))):
                                                  -int(np.ceil(end_skipT*(229390681/4e6)))].quantile(q=0.9)

    # Calculate overall ISO 532-1 loudness from 2-channel time-varing loudness
    loudnessISO1TimeVar = pd.DataFrame(workbookdata['Sheet8'].iloc[13:, 1:3].values,
                                       columns=workbookdata['Sheet8'].iloc[12, 1:3],
                                       index=workbookdata['Sheet8'].iloc[13:, 0])
    # 2-channel overall (5% exceeded = 95th percentile) loudness
    loudnessISO105Ex = loudnessISO1TimeVar.iloc[int(start_skipT*1e3):-int(end_skipT*1e3)].quantile(q=0.95)
    # max of l/r channel (5% exceeded = 95th percentile) loudness
    loudnessISO105ExMaxLR = loudnessISO105Ex.max()

    # 2-channel overall (power-averaged) loudness
    loudnessISO1PowAvg = ((loudnessISO1TimeVar.iloc[int(start_skipT*500):
                                                    -int(end_skipT*500)]**(1/np.log10(2))).sum(axis=0)
                          / len(loudnessISO1TimeVar.iloc[int(start_skipT*500):
                                                         -int(end_skipT*500)]))**np.log10(2)
    # max of l/r channel (95th-percentile) loudness
    loudnessISO1PowAvgMaxLR = loudnessISO1PowAvg.max()

    # Calculate overall ISO 532-3 loudness from binaural time-varing loudness
    loudnessISO3TimeVar = pd.DataFrame(workbookdata['Sheet9'].iloc[13:, 1:2].values,
                                       columns=workbookdata['Sheet9'].iloc[12, 1:2],
                                       index=workbookdata['Sheet9'].iloc[13:, 0])
    # binaural overall (power-averaged) loudness
    loudnessISO3PowAvg = ((loudnessISO3TimeVar.iloc[int(start_skipT*1e3):
                                                    -int(end_skipT*1e3)]**(1/np.log10(2))).sum(axis=0)
                          / len(loudnessISO3TimeVar.iloc[int(start_skipT*1e3):
                                                         -int(end_skipT*1e3)]))**np.log10(2)
    loudnessISO3PowAvg = loudnessISO3PowAvg.iloc[0]  # convert series to float

    # Calculate overall Aures+ISO532-3 sharpness from 2-channel time-varying
    # sharpness
    sharpnessAISO3TimeVar = pd.DataFrame(workbookdata['Sheet10'].iloc[13:, 1:3].values,
                                         columns=workbookdata['Sheet10'].iloc[12, 1:3],
                                         index=workbookdata['Sheet10'].iloc[13:, 0])
    # 2-channel overall (mean) sharpness
    sharpnessAISO3Mean = sharpnessAISO3TimeVar.mean(axis=0)
    # max of l/r channel overall (mean) sharpness
    sharpnessAISO3MeanMaxLR = sharpnessAISO3Mean.max()
    # 2-channel overall (power-averaged) sharpness
    sharpnessAISO3PowAvg = ((sharpnessAISO3TimeVar.iloc[int(start_skipT*1e3):
                                                        -int(end_skipT*1e3)]**(1/np.log10(2))).sum(axis=0)
                            / len(sharpnessAISO3TimeVar.iloc[int(start_skipT*1e3):
                                                             -int(end_skipT*1e3)]))**np.log10(2)
    # max of l/r channel overall (power-averaged) sharpness
    sharpnessAISO3PowAvgMaxLR = sharpnessAISO3PowAvg.max()

    # Calculate overall Sottek Hearing Model impulsiveness from 2-channel
    # time-varying impulsiveness
    impulsivenessHMSTimeVar = pd.DataFrame(workbookdata['Sheet12'].iloc[13:, 1:3].values,
                                           columns=workbookdata['Sheet12'].iloc[12, 1:3],
                                           index=workbookdata['Sheet12'].iloc[13:, 0])
    # 2-channel overall (power-averaged) impulsiveness
    impulsivenessHMSPowAvg = ((impulsivenessHMSTimeVar.iloc[int(np.ceil(60000/55*start_skipT)):
                                                            -int(np.ceil(60000/55*end_skipT))]**(1/np.log10(2))).sum(axis=0)
                              / len(impulsivenessHMSTimeVar.iloc[int(np.ceil(60000/55*start_skipT)):-int(np.ceil(60000/55*end_skipT))]))**np.log10(2)
    # max of l/r overall (power-averaged) impulsiveness
    impulsivenessHMSPowAvgMaxLR = impulsivenessHMSPowAvg.max()
    # 2-channel overall (mean) impulsiveness
    impulsivenessHMSMean = impulsivenessHMSTimeVar.mean(axis=0)
    # max of l/r channel overall (mean) impulsiveness
    impulsivenessHMSMeanMaxLR = impulsivenessHMSMean.max()
    # max of l/r channel overall (max) impulsiveness
    impulsivenessHMSMaxMaxLR = impulsivenessHMSTimeVar.max(axis=None)

    # add results to output DataFrame
    dataByStim.loc[stimulus, 'LoudECMAHMSPowAvgBin'] = loudnessHMSPowAvg
    dataByStim.loc[stimulus, 'TonalECMAHMSAvgMaxLR'] = tonalityHMSAvgMaxLR
    dataByStim.loc[stimulus, 'RoughECMAHMS10ExBin'] = roughnessHMS10Ex
    dataByStim.loc[stimulus, 'FluctstrHMS10ExBin'] = fluctStrHMSTime10Ex
    dataByStim.loc[stimulus, 'LoudISO105ExMaxLR'] = loudnessISO105ExMaxLR
    dataByStim.loc[stimulus, 'LoudISO1PowAvgMaxLR'] = loudnessISO1PowAvgMaxLR
    dataByStim.loc[stimulus, 'LoudISO3PowAvgMaxLR'] = loudnessISO3PowAvg
    dataByStim.loc[stimulus, 'SharpAuresISO3AvgMaxLR'] = sharpnessAISO3MeanMaxLR
    dataByStim.loc[stimulus, 'SharpAuresISO3PowAvgMaxLR'] = sharpnessAISO3PowAvgMaxLR
    dataByStim.loc[stimulus, 'ImpulsHMSPowAvgMaxLR'] = impulsivenessHMSPowAvgMaxLR
    dataByStim.loc[stimulus, 'ImpulsHMSAvgMaxLR'] = impulsivenessHMSMeanMaxLR
    dataByStim.loc[stimulus, 'ImpulsHMSMaxMaxLR'] = impulsivenessHMSMaxMaxLR

# end of for loop over psychoacoustic metrics xlsx files

# add UAS-only and ambient-only data to combined stimuli
# ------------------------------------------------------
indices = indicesAcoustic + indicesPsycho

# Part A UAS only
mask = dataByStim.index.str.find("A_") != -1
UASonlyPtA1 = dataByStim.loc[mask, indices].copy()
UASonlyPtA1 = UASonlyPtA1.add_prefix("UAS", axis=1)
UASonlyPtA2 = UASonlyPtA1.copy()
UASonlyPtA1['newindex'] = [file.replace("A_", "A1_")
                           for file in list(UASonlyPtA1.index)]
UASonlyPtA1.set_index('newindex', inplace=True)
UASonlyPtA2['newindex'] = [file.replace("A_", "A2_")
                           for file in list(UASonlyPtA2.index)]
UASonlyPtA2.set_index('newindex', inplace=True)

UASonlyPtA = pd.concat([UASonlyPtA1, UASonlyPtA2], axis=0)

# Part B UAS only
mask = dataByStim.index.str.find("B_") != -1
UASonlyPtB = dataByStim.loc[mask, indices].copy()
UASonlyPtB = UASonlyPtB.add_prefix("UAS", axis=1)
UASonlyPtB['newindex'] = [file.replace("B_", "B2_")
                          for file in list(UASonlyPtB.index)]
UASonlyPtB.set_index('newindex', inplace=True)

# concatenate UAS only
UASonly = pd.concat([UASonlyPtA, UASonlyPtB], axis=0)

# Part A1 ambient only
mask = dataByStim.index.str.find("A1_CALBIN") != -1
AmbonlyPtA1 = dataByStim.loc[mask, indices].copy()
AmbonlyPtA1 = AmbonlyPtA1.add_prefix("Amb", axis=1)
AmbonlyPtA1 = AmbonlyPtA1.loc[AmbonlyPtA1.index.repeat(
                              sum(dataByStim.index.str.find("A_") != -1) + 1)]
AmbonlyPtA1['newindex'] = UASonlyPtA1.index.union(dataByStim.index[dataByStim.index.str.find("A1_CALBIN") != -1])
AmbonlyPtA1.set_index('newindex', inplace=True)

# Part A2 ambient only
mask = dataByStim.index.str.find("A2_CALBIN") != -1
AmbonlyPtA2 = dataByStim.loc[mask, indices].copy()
AmbonlyPtA2 = AmbonlyPtA2.add_prefix("Amb", axis=1)
AmbonlyPtA2 = AmbonlyPtA2.loc[AmbonlyPtA2.index.repeat(
                              sum(dataByStim.index.str.find("A_") != -1) + 1)]
AmbonlyPtA2['newindex'] = UASonlyPtA2.index.union(dataByStim.index[dataByStim.index.str.find("A2_CALBIN") != -1])
AmbonlyPtA2.set_index('newindex', inplace=True)

AmbonlyPtA = pd.concat([AmbonlyPtA1, AmbonlyPtA2], axis=0)

# Part B ambient only
mask = dataByStim.index.str.find("B2_CALBIN") != -1
AmbonlyPtB = dataByStim.loc[mask, indices].copy()
AmbonlyPtB = AmbonlyPtB.add_prefix("Amb", axis=1)
AmbonlyPtB = AmbonlyPtB.loc[AmbonlyPtB.index.repeat(
                            sum(dataByStim.index.str.find("B_") != -1) + 1)]
AmbonlyPtB['newindex'] = UASonlyPtB.index.union(dataByStim.index[dataByStim.index.str.find("B2_CALBIN") != -1])
AmbonlyPtB.set_index('newindex', inplace=True)

# concatenate UAS only
Ambonly = pd.concat([AmbonlyPtA, AmbonlyPtB], axis=0)

# concatenate UAS only and ambient only
UASAmbonly = pd.concat([UASonly, Ambonly], axis=1)

# merge into output
dataByStim = dataByStim.merge(UASAmbonly.astype(float), how='outer',
                              left_index=True, right_index=True)

# calculate level differences and ratios between UAS and ambient
dataByStim['UASLAeqdiffAmbLAF90'] = dataByStim['UASLAeqMaxLR'] - dataByStim['AmbLAF90ExMaxLR']
dataByStim['UASLASmaxdiffAmbLAF90'] = dataByStim['UASLASmaxMaxLR'] - dataByStim['AmbLAF90ExMaxLR']
dataByStim['UASLASmaxdiffAmbLAF50'] = dataByStim['UASLASmaxMaxLR'] - dataByStim['AmbLAF50ExMaxLR']


# ---------------------
# Partial loudness data
# ---------------------

# open csv file selection dialog and assign filepaths to list
# PROJECT NOTE: the calculated files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Analysis/deeuu_loudness/output?csf=1&web=1&e=ZvblMt
app = QApplication(sys.argv)
fileExts = "*ShortTermPartialLoudness.pkl"
filelist = list(QFileDialog.getOpenFileNames(filter=fileExts))[0]
filelist.sort()

filenames = [filepath.split('/')[-1] for filepath in filelist]

# setup results DataFrame
partLoudness = pd.DataFrame(index=filenames, columns=['UASPartLoudGMSTPowAvg'])

# the deeuu Glasberg&Moore loudness outputs are all sampled at 1 ms
# in each of the deeuu loudness outputs, the first 0.031 seconds is
# redundant 'negative time' and the final 0.169 seconds is redundant added time
# skip indices
skip_starti = 31
skip_endi = -169

# loop over files to extract data
for ii, file in enumerate(filelist):
    
    if filenames[ii].find("A") == 0:
        timeTotal = 25.0
    elif filenames[ii].find("B") == 0:
        timeTotal = 75.0

    partN = pd.read_pickle(file)
    partN = pd.Series(data=partN[skip_starti:skip_endi, 1],
                      index=np.arange(0, timeTotal, 1e-3))

    partialLoudnessPowAvg = (np.sum(partN**(1/np.log10(2)), axis=0)
                             / len(partN))**np.log10(2)

    partLoudness.loc[filenames[ii], 'UASPartLoudGMSTPowAvg'] = partialLoudnessPowAvg


# reindex partial loudness DataFrame to match recording files and merge into
# output DataFrame
partLoudness['newindex'] = [file.replace("ShortTermPartialLoudness", "")
                            for file in list(partLoudness.index)]
partLoudness['newindex'] = [file.replace(".pkl", ".wav")
                            for file in partLoudness['newindex']]

partLoudness.set_index('newindex', inplace=True)

dataByStim = dataByStim.merge(partLoudness.astype(float), how='outer',
                              left_index=True, right_index=True)


# -------------
# Response data
# -------------

# open xlsx file selection dialog and assign filepath
# PROJECT NOTE: the results files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Test_files/Results?csf=1&web=1&e=f3ejcD
app = QApplication(sys.argv)
fileExts = "*.xlsx"
filepath = QFileDialog.getOpenFileName(filter=fileExts)[0]


# Part A
# ------

# read in data and add column indicating the stimulus recording file
partAResponses = pd.read_excel(io=filepath, sheet_name="PtAResponse",
                               header=0)
partAResponses['Recording'] = [StimFile.split('.')[0] + "_CALBIN_Pa.wav"
                               for StimFile in partAResponses['Stim File']]

partAData = partAResponses.drop(columns=['Typology', 'Aircraft', 'Other'])
partAData.columns = [s.replace(' ', '') for s in partAData.columns]

# initialise DataFrames for loop over stimulus recordings
partA = pd.DataFrame(index=partAResponses['Recording'].unique())
partAstats = pd.DataFrame(index=partAResponses['Recording'].unique())

# loop over stimuli recording names, extract corresponding response data and
# tranpose data with recording as index
# apply basic normality tests and calculate aggregate statistics
for ii, file in enumerate(partAResponses['Recording'].unique()):
    partAValence = partAResponses.loc[partAResponses['Recording'] == file,
                                      ['ID#', 'Valence']]
    columns = ["Valence_" + str(ID) for ID in partAValence['ID#']]
    partAValence = pd.DataFrame(data=np.array(partAValence['Valence']),
                                index=columns, columns=[file]).transpose()
    partAArousal = partAResponses.loc[partAResponses['Recording'] == file,
                                      ['ID#', 'Arousal']]
    columns = ["Arousal_" + str(ID) for ID in partAArousal['ID#']]
    partAArousal = pd.DataFrame(data=np.array(partAArousal['Arousal']),
                                index=columns, columns=[file]).transpose()
    partAAnnoy = partAResponses.loc[partAResponses['Recording'] == file,
                                    ['ID#', 'Annoyance']]
    columns = ["Annoyance_" + str(ID) for ID in partAAnnoy['ID#']]
    partAAnnoy = pd.DataFrame(data=np.array(partAAnnoy['Annoyance']),
                              index=columns, columns=[file]).transpose()
    partANotice = partAResponses.loc[partAResponses['Recording'] == file,
                                     ['ID#', 'UAS_noticed']]
    columns = ["UAS_noticed_" + str(ID) for ID in partANotice['ID#']]
    partANotice = pd.DataFrame(data=np.array(partANotice['UAS_noticed']),
                               index=columns, columns=[file]).transpose()       

    # response data normality testing and aggregation
    # (partAValence.transpose()).plot.hist(bins=np.arange(0.5, 7.5, 1),
    #                                      alpha=0.75, xlabel="Valence rating")
    # (partAArousal.transpose()).plot.hist(bins=np.arange(0.5, 7.5, 1),
    #                                      alpha=0.75, xlabel="Arousal rating")
    # (partAAnnoy.transpose()).plot.hist(bins=np.arange(-0.5, 11.5, 1),
    #                                    alpha=0.75, xlabel="Annoyance rating")
    valenceSWtest = pd.DataFrame(data=[np.array(stats.shapiro(partAValence))],
                                 columns=['ValenceSWtestW', 'ValenceSWtestp'],
                                 index=[file])
    arousalSWtest = pd.DataFrame(data=[np.array(stats.shapiro(partAArousal))],
                                 columns=['ArousalSWtestW', 'ArousalSWtestp'],
                                 index=[file])
    annoySWtest = pd.DataFrame(data=[np.array(stats.shapiro(partAAnnoy))],
                               columns=['AnnoySWtestW', 'AnnoySWtestp'],
                               index=[file])
    valenceAgg = pd.DataFrame(data=[[np.percentile(partAValence.values,
                                                   q=50, axis=1)[0],
                                    np.mean(partAValence.values, axis=1)[0]]],
                              columns=['ValenceMedian', 'ValenceMean'],
                              index=[file])
    arousalAgg = pd.DataFrame(data=[[np.percentile(partAArousal.values,
                                                   q=50, axis=1)[0],
                                    np.mean(partAArousal.values, axis=1)[0]]],
                              columns=['ArousalMedian', 'ArousalMean'],
                              index=[file])
    annoyAgg = pd.DataFrame(data=[[np.percentile(partAAnnoy.values,
                                                 q=50, axis=1)[0],
                                  np.mean(partAAnnoy.values, axis=1)[0]]],
                            columns=['AnnoyMedian', 'AnnoyMean'],
                            index=[file])
    noticeAgg = pd.DataFrame(data=[[np.sum(partANotice.values, axis=1)[0],
                                   np.mean(partANotice.values, axis=1)[0]]],
                             columns=['NoticedTotal', 'NoticedProportion'],
                             index=[file])

    # add results to DataFrame
    if ii == 0:
        partA = partA.join([partAValence, partAArousal, partAAnnoy,
                            partANotice, valenceAgg, arousalAgg,
                            annoyAgg, noticeAgg])
        partAstats = partAstats.join([valenceSWtest, arousalSWtest,
                                      annoySWtest, valenceAgg, arousalAgg,
                                      annoyAgg, noticeAgg])
    else:
        partAstats.loc[file, valenceSWtest.columns] = valenceSWtest.loc[file]
        partAstats.loc[file, arousalSWtest.columns] = arousalSWtest.loc[file]
        partAstats.loc[file, annoySWtest.columns] = annoySWtest.loc[file]
        partAstats.loc[file, valenceAgg.columns] = valenceAgg.loc[file]
        partAstats.loc[file, arousalAgg.columns] = arousalAgg.loc[file]
        partAstats.loc[file, annoyAgg.columns] = annoyAgg.loc[file]
        partAstats.loc[file, noticeAgg.columns] = noticeAgg.loc[file]
    
        partA.loc[file, partAValence.columns] = partAValence.loc[file]
        partA.loc[file, partAArousal.columns] = partAArousal.loc[file]
        partA.loc[file, partAAnnoy.columns] = partAAnnoy.loc[file]
        partA.loc[file, partANotice.columns] = partANotice.loc[file]
        partA.loc[file, valenceAgg.columns] = valenceAgg.loc[file]
        partA.loc[file, arousalAgg.columns] = arousalAgg.loc[file]
        partA.loc[file, annoyAgg.columns] = annoyAgg.loc[file]
        partA.loc[file, noticeAgg.columns] = noticeAgg.loc[file]


# Part B
# ------

# read in data and add column indicating the stimulus recording file
partBResponses = pd.read_excel(io=filepath, sheet_name="PtBResponse",
                               header=0)
partBResponses['Recording'] = [StimFile.split('.')[0] + "_CALBIN_Pa.wav"
                                for StimFile in partBResponses['Stim File']]

partBData = partBResponses.copy()
partBData.columns = [s.replace(' ', '') for s in partBData.columns]

# initialise DataFrames for loop over stimulus recordings
partB = pd.DataFrame(index=partBResponses['Recording'].unique())
partBstats = pd.DataFrame(index=partBResponses['Recording'].unique())

# loop over stimuli recording names, extract corresponding response data and
# tranpose data with recording as index
# apply basic normality tests and calculate aggregate statistics
for ii, file in enumerate(partBResponses['Recording'].unique()):
    partBValence = partBResponses.loc[partBResponses['Recording'] == file,
                                      ['ID#', 'Valence']]
    columns = ["Valence_" + str(ID) for ID in partBValence['ID#']]
    partBValence = pd.DataFrame(data=np.array(partBValence['Valence']),
                                index=columns, columns=[file]).transpose()
    partBArousal = partBResponses.loc[partBResponses['Recording'] == file,
                                      ['ID#', 'Arousal']]
    columns = ["Arousal_" + str(ID) for ID in partBArousal['ID#']]
    partBArousal = pd.DataFrame(data=np.array(partBArousal['Arousal']),
                                index=columns, columns=[file]).transpose()
    partBAnnoy = partBResponses.loc[partBResponses['Recording'] == file,
                                    ['ID#', 'Annoyance']]
    columns = ["Annoyance_" + str(ID) for ID in partBAnnoy['ID#']]
    partBAnnoy = pd.DataFrame(data=np.array(partBAnnoy['Annoyance']),
                              index=columns, columns=[file]).transpose()

    # response data normality testing and aggregation
    # (partBValence.transpose()).plot.hist(bins=np.arange(0.5, 7.5, 1),
    #                                      alpha=0.75, xlabel="Valence rating")
    # (partBArousal.transpose()).plot.hist(bins=np.arange(0.5, 7.5, 1),
    #                                      alpha=0.75, xlabel="Arousal rating")
    # (partBAnnoy.transpose()).plot.hist(bins=np.arange(-0.5, 11.5, 1),
    #                                    alpha=0.75, xlabel="Annoyance rating")
    valenceSWtest = pd.DataFrame(data=[np.array(stats.shapiro(partBValence))],
                                 columns=['ValenceSWtestW', 'ValenceSWtestp'],
                                 index=[file])
    arousalSWtest = pd.DataFrame(data=[np.array(stats.shapiro(partBArousal))],
                                 columns=['ArousalSWtestW', 'ArousalSWtestp'],
                                 index=[file])
    annoySWtest = pd.DataFrame(data=[np.array(stats.shapiro(partBAnnoy))],
                               columns=['AnnoySWtestW', 'AnnoySWtestp'],
                               index=[file])
    valenceAgg = pd.DataFrame(data=[[np.percentile(partBValence.values,
                                                   q=50, axis=1)[0],
                                    np.mean(partBValence.values, axis=1)[0]]],
                              columns=['ValenceMedian', 'ValenceMean'],
                              index=[file])
    arousalAgg = pd.DataFrame(data=[[np.percentile(partBArousal.values,
                                                   q=50, axis=1)[0],
                                    np.mean(partBArousal.values, axis=1)[0]]],
                              columns=['ArousalMedian', 'ArousalMean'],
                              index=[file])
    annoyAgg = pd.DataFrame(data=[[np.percentile(partBAnnoy.values,
                                                 q=50, axis=1)[0],
                                  np.mean(partBAnnoy.values, axis=1)[0]]],
                            columns=['AnnoyMedian', 'AnnoyMean'],
                            index=[file])

    # add results to DataFrame
    if ii == 0:
        partB = partB.join([partBValence, partBArousal, partBAnnoy,
                            valenceAgg, arousalAgg, annoyAgg])
        partBstats = partBstats.join([valenceSWtest, arousalSWtest,
                                      annoySWtest, valenceAgg, arousalAgg,
                                      annoyAgg])
    else:
        partBstats.loc[file, valenceSWtest.columns] = valenceSWtest.loc[file]
        partBstats.loc[file, arousalSWtest.columns] = arousalSWtest.loc[file]
        partBstats.loc[file, annoySWtest.columns] = annoySWtest.loc[file]
        partBstats.loc[file, valenceAgg.columns] = valenceAgg.loc[file]
        partBstats.loc[file, arousalAgg.columns] = arousalAgg.loc[file]
        partBstats.loc[file, annoyAgg.columns] = annoyAgg.loc[file]
        
        partB.loc[file, partBValence.columns] = partBValence.loc[file]
        partB.loc[file, partBArousal.columns] = partBArousal.loc[file]
        partB.loc[file, partBAnnoy.columns] = partBAnnoy.loc[file]
        partB.loc[file, valenceAgg.columns] = valenceAgg.loc[file]
        partB.loc[file, arousalAgg.columns] = arousalAgg.loc[file]
        partB.loc[file, annoyAgg.columns] = annoyAgg.loc[file]

allResponses = pd.concat([partB, partA], axis=0, join='outer')
allResponses.sort_index(inplace=True)

# merge response data into output
dataByStim = dataByStim.merge(allResponses, how='outer',
                              left_index=True, right_index=True)

# add pre- and post-test response responses to dataset
preTestResponses = pd.read_excel(io=filepath, sheet_name="Pre-test",
                                 header=0)
preTestResponses.drop(columns=['Impairment_details'], inplace=True)
preTestResponses.drop(columns=preTestResponses.loc[:,
                                                   'PANAS_interested':].columns,
                      inplace=True)
postTestResponses = pd.read_excel(io=filepath, sheet_name="Post-test",
                                  header=0)
postTestResponses.drop(columns=['Nationality', 'Language',
                                'Response Considerations',
                                'Source comments', 'AAM Exp',
                                'General Feedback'], inplace=True)
postTestResponses.drop(columns=postTestResponses.loc[:, 'NSS21':'NSS8'].columns,
                       inplace=True)

prePostTestResponses = preTestResponses.merge(postTestResponses, on='ID#')

# -------------------------------------
# Add categorical variables for stimuli
# -------------------------------------

# NOTE: the order of this section matters, as some categorical variable
# definitions depend on Boolean logic derived from other categorical variables
# added

# Experiment session part A or B
dataByStim.loc[np.logical_or(np.logical_or(dataByStim.index.str.find("A1_")
                                           == 0,
                                           dataByStim.index.str.find("A2_")
                                           == 0),
                             dataByStim.index.str.find("A_") == 0),
               'SessionPart'] = "A"
dataByStim.loc[np.logical_or(np.logical_or(dataByStim.index.str.find('B1_')
                                           == 0,
                                           dataByStim.index.str.find('B2_')
                                           == 0),
                             dataByStim.index.str.find("B_") == 0),
               'SessionPart'] = "B"

dataByStim.insert(loc=0, column='SessionPart',
                  value=dataByStim.pop('SessionPart'))

# Stimulus duration
for file in dataByStim.index:
    if file.split('_')[0].find("A") == 0:
        dataByStim.loc[file, 'StimDuration'] = 25
    elif file.split('_')[0].find("B") == 0:
        dataByStim.loc[file, 'StimDuration'] = 75
    else:
        dataByStim.loc[file, 'StimDuration'] = np.nan

dataByStim.insert(loc=1, column='StimDuration',
                  value=dataByStim.pop('StimDuration'))

# UAS LAeq level category
# Add the No UAS entries first
dataByStim.loc[np.logical_or(dataByStim.index.str.find("B2_CALBIN") != -1,
                             np.logical_or(dataByStim.index.str.find("A1_CALBIN") != -1,
                                           dataByStim.index.str.find("A2_CALBIN") != -1)),
               'UASLAeq'] = "No UAS"

dataByStim.loc[np.logical_or(dataByStim.index.str.find("_F_1") != -1,
                             np.logical_and(dataByStim.index.str.find("_1_")
                                            != -1,
                                            dataByStim.index.str.find("_F_2")
                                            == -1)),
               'UASLAeq'] = int(60)
dataByStim.loc[np.logical_or(dataByStim.index.str.find("_F_2") != -1,
                             np.logical_and(dataByStim.index.str.find("_2_")
                                            != -1,
                                            dataByStim.index.str.find("_F_1")
                                            == -1)),
               'UASLAeq'] = int(54)
dataByStim.loc[np.logical_and(dataByStim.index.str.find("_3_") != -1,
                              ~np.logical_or(dataByStim.index.str.find("_F_1")
                                             != -1,
                                             dataByStim.index.str.find("_F_2")
                                             != -1)),
               'UASLAeq'] = int(48)
dataByStim.loc[dataByStim.index.str.find("_4_") != -1, 'UASLAeq'] = int(42)                           

dataByStim.insert(loc=2, column='UASLAeq', value=dataByStim.pop('UASLAeq'))

# ambient sound LAeq level category
dataByStim.loc[dataByStim.index.str.find("A1") != -1, 'AmbientLAeq'] = int(58)
dataByStim.loc[dataByStim.index.str.find("A2") != -1, 'AmbientLAeq'] = int(52)
dataByStim.loc[np.logical_or(dataByStim.index.str.find("B1") != -1,
                             dataByStim.index.str.find("B2") != -1),
               'AmbientLAeq'] = int(52)
dataByStim.insert(loc=3, column='AmbientLAeq',
                  value=dataByStim.pop('AmbientLAeq'))

# UAS SNR category
dataByStim.loc[np.logical_or(dataByStim.index.str.find("B2_CALBIN") != -1,
                             np.logical_or(dataByStim.index.str.find("A1_CALBIN") != -1,
                                           dataByStim.index.str.find("A2_CALBIN") != -1)),
               'SNRlevel'] = "No UAS"
mask = np.logical_and(np.logical_or
                      (np.logical_or(dataByStim.index.str.find("A1") != -1,
                                     dataByStim.index.str.find("A2") != -1),
                       np.logical_or(dataByStim.index.str.find("B1") != -1,
                                     dataByStim.index.str.find("B2") != -1)),
                      ~np.logical_or
                      (np.logical_or
                       (dataByStim.index.str.find("A1_CALBIN") != -1,
                        dataByStim.index.str.find("A2_CALBIN") != -1),
                       dataByStim.index.str.find("B2_CALBIN") != -1))
dataByStim.loc[mask, 'SNRlevel'] = (dataByStim.loc[mask, 'UASLAeq']
                                    - dataByStim.loc[mask, 'AmbientLAeq']).astype(int)
dataByStim.insert(loc=4, column='SNRlevel', value=dataByStim.pop('SNRlevel'))

# recast UASLAeq and SNRlevel as ordered categorical types
dataByStim['SNRlevel'] = dataByStim['SNRlevel'].astype(str)
dataByStim['UASLAeq'] = dataByStim['UASLAeq'].astype(str)
dataByStim['SNRlevel'] = pd.Categorical(dataByStim['SNRlevel'], ["No UAS",
                                                                 "-16", "-10",
                                                                 "-4", "2",
                                                                 "8"])
dataByStim['UASLAeq'] = pd.Categorical(dataByStim['UASLAeq'], ["No UAS", "42", "48",
                                                               "54", "60"])

# UAS operation
dataByStim.loc[dataByStim.index.str.find("_F_") != -1, 'UASOperation'] = "Flyby"
dataByStim.loc[dataByStim.index.str.find("_T_") != -1, 'UASOperation'] = "Takeoff"
dataByStim.loc[np.logical_and(dataByStim.index.str.find("_L_") != -1,
                              dataByStim.index.str.find("_F_") == -1),
                              'UASOperation'] = "Landing"
dataByStim.loc[np.logical_or(dataByStim.index.str.find("B2_CALBIN") != -1,
                             np.logical_or(dataByStim.index.str.find("A1_CALBIN")
                                           != -1,
                                           dataByStim.index.str.find("A2_CALBIN")
                                           != -1)),
               'UASOperation'] = "No UAS"

dataByStim.insert(loc=5, column='UASOperation',
                  value=dataByStim.pop('UASOperation'))

# UAS event quantity / density
dataByStim.loc[dataByStim['UASOperation'] == "No UAS", 'UASEvents'] = 0
dataByStim.loc[np.logical_and(dataByStim['SessionPart'] == "A",
                              dataByStim['UASOperation'] != "No UAS"),
               'UASEvents'] = 1
dataByStim.loc[np.logical_and(dataByStim['SessionPart'] == "B",
                              dataByStim.index.str.find("_1_CALBIN") != -1),
               'UASEvents'] = 1
dataByStim.loc[np.logical_and(dataByStim['SessionPart'] == "B",
                              dataByStim.index.str.find("_3_CALBIN") != -1),
               'UASEvents'] = 3
dataByStim.loc[np.logical_and(dataByStim['SessionPart'] == "B",
                              dataByStim.index.str.find("_5_CALBIN") != -1),
               'UASEvents'] = 5
dataByStim.loc[np.logical_and(dataByStim['SessionPart'] == "B",
                              dataByStim.index.str.find("_9_CALBIN") != -1),
               'UASEvents'] = 9

dataByStim.insert(loc=6, column='UASEvents', value=dataByStim.pop('UASEvents'))

dataByStim['UASEventDensity'] = dataByStim['UASEvents']/dataByStim['StimDuration']

dataByStim.insert(loc=7, column='UASEventDensity',
                  value=dataByStim.pop('UASEventDensity'))

# UAS type
dataByStim.loc[dataByStim.index.str.find("_YnTy_") != -1, 'UASType'] = "YnTy"
dataByStim.loc[dataByStim.index.str.find("_M300_") != -1, 'UASType'] = "M300"
dataByStim.loc[dataByStim.index.str.find("_T150_") != -1, 'UASType'] = "T150"
dataByStim.loc[np.logical_or(dataByStim.index.str.find("B2_CALBIN") != -1,
                             np.logical_or(dataByStim.index.str.find("A1_CALBIN") != -1,
                                           dataByStim.index.str.find("A2_CALBIN") != -1)),
               'UASType'] = "No UAS"
dataByStim['UASType'] = pd.Categorical(dataByStim['UASType'], ["No UAS",
                                                               "YnTy",
                                                               "M300",
                                                               "T150"])
dataByStim.insert(loc=8, column='UASType',
                  value=dataByStim.pop('UASType'))

# ambient sound environment
dataByStim.loc[np.logical_or(dataByStim.index.str.find("A2") != -1,
                             dataByStim.index.str.find("B2") != -1),
               'AmbientEnv'] = "Park"
dataByStim.loc[dataByStim.index.str.find("A1") != -1, 'AmbientEnv'] = "Street"

dataByStim.insert(loc=9, column='AmbientEnv',
                  value=dataByStim.pop('AmbientEnv'))


# separate 'by stimulus' output into test data and auxiliary data and save to
# file
dataByStimTest = dataByStim.loc[np.logical_or(
                                dataByStim.index.str.find("B2") == 0,
                                np.logical_or(dataByStim.index.str.find("A1")
                                              == 0,
                                              dataByStim.index.str.find("A2")
                                              == 0)), :]
dataByStimTestA = dataByStimTest.loc[np.logical_or(
                                     dataByStimTest.index.str.find("A1") == 0,
                                     dataByStimTest.index.str.find("A2") == 0),
                                     :].dropna(axis=1, how='all')
dataByStimTestB = dataByStimTest.loc[dataByStimTest.index.str.find("B2") == 0,
                                     :].dropna(axis=1, how='all')

dataByStimAux = dataByStim.loc[~np.logical_or(
                               dataByStim.index.str.find("B2") == 0,
                               np.logical_or(dataByStim.index.str.find("A1")
                                             == 0,
                                             dataByStim.index.str.find("A2")
                                             == 0)),
                               :].copy().dropna(axis=1, how='all')

app = QApplication(sys.argv)
outFilePath = QFileDialog.getExistingDirectory()

dataByStimTestFilePath = os.path.join(outFilePath,
                                      "refmap_listest1_testdata_ByStim.csv")
dataByStimTest.to_csv(dataByStimTestFilePath)
dataByStimTestAFilePath = os.path.join(outFilePath,
                                       "refmap_listest1_testdataA_ByStim.csv")
dataByStimTestA.to_csv(dataByStimTestAFilePath)
dataByStimTestBFilePath = os.path.join(outFilePath,
                                       "refmap_listest1_testdataB_ByStim.csv")
dataByStimTestB.to_csv(dataByStimTestBFilePath)
dataByStimAuxFilePath = os.path.join(outFilePath,
                                     "refmap_listest1_auxdata.csv")
dataByStimAux.to_csv(dataByStimAuxFilePath)
dataByStimFilePath = os.path.join(outFilePath,
                                  "refmap_listest1_alldata_ByStim.csv")
dataByStim.to_csv(dataByStimFilePath)

# form sub datasets for inter-rater reliability testing
subdataByStimAnnoy = dataByStim.loc[:, "Annoyance_1":"Annoyance_19"].copy()
subdataByStimValence = dataByStim.loc[:, "Valence_1":"Valence_19"].copy()
subdataByStimArousal = dataByStim.loc[:, "Arousal_1":"Arousal_19"].copy()

subdataByStimAnnoy.set_index(dataByStim.index, inplace=True)
subdataByStimAnnoy = subdataByStimAnnoy.transpose()
subdataByStimAnnoyFilePath = os.path.join(outFilePath,
                                          "refmap_listest1_testdataAnnoyRater.csv")
subdataByStimAnnoy.to_csv(subdataByStimAnnoyFilePath)

subdataByStimValence.set_index(dataByStim.index, inplace=True)
subdataByStimValence = subdataByStimValence.transpose()
subdataByStimValenceFilePath = os.path.join(outFilePath,
                                          "refmap_listest1_testdataValenceRater.csv")
subdataByStimValence.to_csv(subdataByStimValenceFilePath)

subdataByStimArousal.set_index(dataByStim.index, inplace=True)
subdataByStimArousal = subdataByStimArousal.transpose()
subdataByStimArousalFilePath = os.path.join(outFilePath,
                                          "refmap_listest1_testdataArousalRater.csv")
subdataByStimArousal.to_csv(subdataByStimArousalFilePath)

# merge response and stimuli data into 'by participant' test datasets, and
# save to file

partADataBySubj = pd.merge(left=partAData,
                           right=dataByStimTestA.loc[:,
                                                     :'UASPartLoudGMSTPowAvg'],
                           how='outer', left_on='Recording', right_index=True)
partADataBySubj.sort_values(by='ID#', axis=0, inplace=True)
partADataBySubj = pd.merge(left=partADataBySubj,
                           right=prePostTestResponses, how='left',
                           left_on='ID#', right_on='ID#')
partADataBySubj.drop(columns='Part', inplace=True)
partADataBySubj.insert(loc=0, column='SessionPart',
                       value=partADataBySubj.pop('SessionPart'))

partBDataBySubj = pd.merge(left=partBData,
                           right=dataByStimTestB.loc[:,
                                                     :'UASPartLoudGMSTPowAvg'],
                           how='left', left_on='Recording', right_index=True)
partBDataBySubj.sort_values(by='ID#', axis=0, inplace=True)
partBDataBySubj = pd.merge(left=partBDataBySubj,
                           right=prePostTestResponses, how='outer',
                           left_on='ID#', right_on='ID#')
partBDataBySubj.drop(columns='Part', inplace=True)
partBDataBySubj.insert(loc=0, column='SessionPart',
                       value=partBDataBySubj.pop('SessionPart'))

allDataBySubj = pd.concat([partADataBySubj, partBDataBySubj], axis=0,
                          ignore_index=True)

allDataFilePath = os.path.join(outFilePath,
                               "refmap_listest1_testdata_BySubj.csv")
allDataBySubj.to_csv(allDataFilePath, index=False)
partADataFilePath = os.path.join(outFilePath,
                                 "refmap_listest1_testdataA_BySubj.csv")
partADataBySubj.to_csv(partADataFilePath, index=False)
partBDataFilePath = os.path.join(outFilePath,
                                 "refmap_listest1_testdataB_BySubj.csv")
partBDataBySubj.to_csv(partBDataFilePath, index=False)

# save pre-test and post-test response data to separate files

preTestDataFilePath = os.path.join(outFilePath,
                                   "refmap_listest1_pretestdata.csv")
preTestResponses.to_csv(preTestDataFilePath, index=False)

postTestDataFilePath = os.path.join(outFilePath,
                                   "refmap_listest1_posttestdata.csv")
postTestResponses.to_csv(postTestDataFilePath, index=False)