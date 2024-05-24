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
import dsp.noct
from scipy import stats, optimize
# import matlab.engine

# create MATLAB engine object
# eng = matlab.engine.start_matlab()

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
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Calibration/Post_calib_recs/Pa_calib?csf=1&web=1&e=7eX89P
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

# -----------------
# EPNL calculations
# -----------------

# # output variables
indicesPNL = ["EPNLMaxLR", "EPNLNoTMaxLR", "PNLTmaxMaxLR"]

# # Select Audio2EPNL function file folder
# Audio2EPNLPath = QFileDialog.getExistingDirectory()
# eng.addpath(eng.genpath(Audio2EPNLPath));
# PNLout= eng.Audio2EPNL(1, filelist, nargout=3)
#TODO
app = QApplication(sys.argv)
fileExts = "*.csv"
filelist = list(QFileDialog.getOpenFileNames(filter=fileExts))[0]
filelist.sort()
EPNL = pd.read_csv(filelist[0], header=0, index_col=0)
EPNL_No_Tone = pd.read_csv(filelist[1], header=0, index_col=0)
PNLTMax = pd.read_csv(filelist[2], header=0, index_col=0)

dataByStim = dataByStim.merge(pd.DataFrame(EPNL.max(axis=1),
                                           columns=[indicesPNL[0]]),
                              how='outer', left_index=True, right_index=True)
dataByStim = dataByStim.merge(pd.DataFrame(EPNL_No_Tone.max(axis=1),
                                           columns=[indicesPNL[1]]),
                              how='outer', left_index=True, right_index=True)
dataByStim = dataByStim.merge(pd.DataFrame(PNLTMax.max(axis=1),
                                           columns=[indicesPNL[2]]),
                              how='outer', left_index=True, right_index=True)

# ------------------------------------------------
# Psychoacoustic metrics single values calculation
# ------------------------------------------------

# output variables
indicesPsycho = ["LoudECMAHMSPowAvgBin",
                  "TonalECMAHMSAvgMaxLR",
                  "TonalECMAHMS05ExMaxLR",
                  "RoughECMAHMS10ExBin",
                  "RoughECMAHMS05ExBin",
                  "FluctstrHMS10ExBin",
                  "FluctstrHMS05ExBin",
                  "LoudISO105ExMaxLR",
                  "LoudISO1PowAvgMaxLR",
                  "LoudISO3PowAvgMaxLR",
                  "SharpAuresISO3AvgMaxLR",
                  "SharpAuresISO3PowAvgMaxLR",
                  "SharpAuresISO305ExMaxLR",
                  "ImpulsHMSPowAvgMaxLR",
                  "ImpulsHMSAvgMaxLR",
                  "ImpulsHMSMaxMaxLR",
                  "ImpulsHMS05ExMaxLR",
                  "dTonalECMAHMSAvgMaxLR",
                  "dTonalECMAHMS05ExMaxLR",
                  "dRoughECMAHMS10ExBin",
                  "dRoughECMAHMS05ExBin",
                  "dFluctstrHMS10ExBin",
                  "dFluctstrHMS05ExBin",
                  "dSharpAuresISO3AvgMaxLR",
                  "dSharpAuresISO3PowAvgMaxLR",
                  "dSharpAuresISO305ExMaxLR",
                  "dImpulsHMSPowAvgMaxLR",
                  "dImpulsHMSAvgMaxLR",
                  "dImpulsHMSMaxMaxLR",
                  "dImpulsHMS05ExMaxLR"]

dataByStim = pd.concat([dataByStim, pd.DataFrame(index=dataByStim.index,
                                                 columns=indicesPsycho,
                                                 dtype=float)], axis=1)

# output SQM sample rates
sampleRateLoudECMA = 187.5
sampleRateTonalECMA = 187.5
sampleRateLoudISO1 = 500
sampleRateLoudISO3 = 1e3
sampleRateRoughECMA = 50
sampleRateFluctHMS = 229390681/4e6
sampleRateImpulsHMS = 60000/55

# time value for moving averaging of SQMs for difference calculations
windowT = 0.05

# open xlsx file selection dialog and assign filepaths to list
# PROJECT NOTE: the calculated files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Analysis/ArtemiS/Output?csf=1&web=1&e=0TAPVz
app = QApplication(sys.argv)
fileExts = "*.xlsx"
filelist = list(QFileDialog.getOpenFileNames(filter=fileExts,
                                             caption=r"Select ArtemiS output files"))[0]
filelist.sort()
filenames = [filepath.split('/')[-1] for filepath in filelist]

# loop over files to analyse
for file in filelist:
    print(file.split('/')[-1])
    workbookdata = pd.read_excel(io=file, sheet_name=None)
    stimulus = workbookdata['Sheet1'].columns[5].split(sep='\'')[1]
    # Calculate ECMA-418-2:2022 Sottek Hearing Model overall loudness from
    # 2-channel specific loudness
    # left channel
    specLoudnessHMSL = pd.DataFrame(workbookdata['Sheet1'].iloc[14:, 1:].values,
                                    columns=workbookdata['Sheet1'].iloc[13, 1:],
                                    index=workbookdata['Sheet1'].iloc[14:, 0])
    # # right channel
    specLoudnessHMSR = pd.DataFrame(workbookdata['Sheet2'].iloc[14:, 1:].values,
                                    columns=workbookdata['Sheet2'].iloc[13, 1:],
                                    index=workbookdata['Sheet2'].iloc[14:, 0])
    # # binaural specific loudness (ECMA-418-2:2022 Equation 118)
    specLoudnessHMSBin = ((specLoudnessHMSL**2
                            + specLoudnessHMSR**2)/2).pow(1/2)
    # # binaural time-dependent loudness (ECMA-418-2:2022 Equation 116)
    loudnessHMSTimeVar = specLoudnessHMSBin.sum(axis=0)*0.5
    # # binaural overall (power-averaged) loudness (ECMA-418-2:2022 Equation 117)
    loudnessHMSPowAvg = ((loudnessHMSTimeVar.iloc[int(np.ceil(sampleRateLoudECMA*start_skipT)):
                                                  -int(np.ceil(sampleRateLoudECMA*end_skipT))]**(1/np.log10(2))).sum()
                          / len(loudnessHMSTimeVar.iloc[int(np.ceil(sampleRateLoudECMA*start_skipT)):
                                                        -int(np.ceil(sampleRateLoudECMA*end_skipT))]))**np.log10(2)

    # Calculate ECMA-418-2:2022 Sottek Hearing Model overall tonality from
    # 2-channel specific tonality
    # left channel
    specTonalityHMSL = pd.DataFrame(workbookdata['Sheet3'].iloc[14:, 1:].values,
                                    columns=workbookdata['Sheet3'].iloc[13, 1:],
                                    index=workbookdata['Sheet3'].iloc[14:, 0])
    # right channel
    specTonalityHMSR = pd.DataFrame(workbookdata['Sheet4'].iloc[14:, 1:].values,
                                    columns=workbookdata['Sheet4'].iloc[13, 1:],
                                    index=workbookdata['Sheet4'].iloc[14:, 0])
    # 2-channel time-varing tonality (max, not integration)
    tonalityHMSTimeVar = pd.concat([specTonalityHMSL.max(axis=0),
                                    specTonalityHMSR.max(axis=0)],
                                    axis=1)
    # 2-channel time-averaged tonality (omitting T<=0.02)
    tonalityHMSAvgL = tonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                              -int(np.ceil(sampleRateTonalECMA*end_skipT)),
                                              0][tonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                                                          -int(np.ceil(sampleRateTonalECMA*end_skipT)), 0]
                                                                          > 0.02].mean(axis=0)
    tonalityHMSAvgR = tonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                              -int(np.ceil(sampleRateTonalECMA*end_skipT)),
                                              1][tonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                                                          -int(np.ceil(sampleRateTonalECMA*end_skipT)), 1]
                                                                          > 0.02].mean(axis=0)
    tonalityHMSAvgMaxLR = max(tonalityHMSAvgL, tonalityHMSAvgR)
    # 2-channel 5% exceeded tonality (omitting T<=0.02)
    tonalityHMS05ExL = tonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                              -int(np.ceil(sampleRateTonalECMA*end_skipT)),
                                              0][tonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                                                          -int(np.ceil(sampleRateTonalECMA*end_skipT)), 0]
                                                                          > 0.02].quantile(q=0.95)
    tonalityHMS05ExR = tonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                                -int(np.ceil(sampleRateTonalECMA*end_skipT)),
                                                1][tonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                                                          -int(np.ceil(sampleRateTonalECMA*end_skipT)), 1]
                                                                          > 0.02].quantile(q=0.95)
    tonalityHMS05ExMaxLR = max(tonalityHMS05ExL, tonalityHMS05ExR)

    # Calculate ECMA-418-2:2022 Sottek Hearing Model binaural overall roughness
    # from binaural specific roughness
    specRoughHMSBin = pd.DataFrame(workbookdata['Sheet5'].iloc[14:, 1:].values,
                                    columns=workbookdata['Sheet5'].iloc[13, 1:],
                                    index=workbookdata['Sheet5'].iloc[14:, 0])
    # binaural time-varying roughness
    roughHMSTimeVar = specRoughHMSBin.sum(axis=0)*0.5
    # binaural overall (90th percentile = 10% exceeded) roughness
    roughHMS10Ex = roughHMSTimeVar.iloc[int(np.ceil(sampleRateRoughECMA*start_skipT)):
                                        -int(np.ceil(sampleRateRoughECMA*end_skipT))].quantile(q=0.9)
    # binaural overall (95th percentile = 5% exceeded) roughness
    roughHMS05Ex = roughHMSTimeVar.iloc[int(np.ceil(sampleRateRoughECMA*start_skipT)):
                                        -int(np.ceil(sampleRateRoughECMA*end_skipT))].quantile(q=0.95)

    # Calculate overall Sottek Hearing Model fluctuation strength from
    # 2-channel specific fluctuation strength 
    # calculate according to ECMA-418-2:2022 approach for roughness
    # left channel
    specFluctStrHMSL = pd.DataFrame(workbookdata['Sheet6'].iloc[14:, 1:].values,
                                    columns=workbookdata['Sheet6'].iloc[13, 1:],
                                    index=workbookdata['Sheet6'].iloc[14:, 0])
    # right channel
    specFluctStrHMSR = pd.DataFrame(workbookdata['Sheet7'].iloc[14:, 1:].values,
                                    columns=workbookdata['Sheet7'].iloc[13, 1:],
                                    index=workbookdata['Sheet7'].iloc[14:, 0])
    # binaural specific fluctuation strength
    # (using ECMA-418-2:2022 Equation 112 for roughness)
    specFluctStrHMSBin = ((specFluctStrHMSL**2
                            + specFluctStrHMSR**2)/2).pow(1/2)
    # binaural time-dependent fluctuation strength
    # (using ECMA-418-2:2022 Equation 111 for roughness)
    fluctStrHMSTimeVar = specFluctStrHMSBin.sum(axis=0)*0.5
    # binaural overall (90th percentile = 10% exceeded) fluctuation strength
    # (using ECMA-418-2:2022 Section 7.1.8 for roughness)
    fluctStrHMSTime10Ex = fluctStrHMSTimeVar.iloc[int(np.ceil(start_skipT*(sampleRateFluctHMS))):
                                                  -int(np.ceil(end_skipT*(sampleRateFluctHMS)))].quantile(q=0.9)
    # binaural overall (95th percentile = 5% exceeded) fluctuation strength
    fluctStrHMSTime05Ex = fluctStrHMSTimeVar.iloc[int(np.ceil(start_skipT*(sampleRateFluctHMS))):
                                                  -int(np.ceil(end_skipT*(sampleRateFluctHMS)))].quantile(q=0.95)
    # # Calculate overall ISO 532-1 loudness from 2-channel time-varing loudness
    loudnessISO1TimeVar = pd.DataFrame(workbookdata['Sheet8'].iloc[13:, 1:3].values,
                                        columns=workbookdata['Sheet8'].iloc[12, 1:3],
                                        index=workbookdata['Sheet8'].iloc[13:, 0])
    # # 2-channel overall (5% exceeded = 95th percentile) loudness
    loudnessISO105Ex = loudnessISO1TimeVar.iloc[int(start_skipT*sampleRateLoudISO1):
                                                -int(end_skipT*sampleRateLoudISO1)].quantile(q=0.95)
    # # max of l/r channel (5% exceeded = 95th percentile) loudness
    loudnessISO105ExMaxLR = loudnessISO105Ex.max()

    # # 2-channel overall (power-averaged) loudness
    loudnessISO1PowAvg = ((loudnessISO1TimeVar.iloc[int(start_skipT*sampleRateLoudISO1):
                                                    -int(end_skipT*sampleRateLoudISO1)]**(1/np.log10(2))).sum(axis=0)
                          / len(loudnessISO1TimeVar.iloc[int(start_skipT*sampleRateLoudISO1):
                                                          -int(end_skipT*sampleRateLoudISO1)]))**np.log10(2)
    # # max of l/r channel (95th-percentile) loudness
    loudnessISO1PowAvgMaxLR = loudnessISO1PowAvg.max()

    # # Calculate overall ISO 532-3 loudness from binaural time-varing loudness
    loudnessISO3TimeVar = pd.DataFrame(workbookdata['Sheet9'].iloc[13:, 1:2].values,
                                        columns=workbookdata['Sheet9'].iloc[12, 1:2],
                                        index=workbookdata['Sheet9'].iloc[13:, 0])
    # # binaural overall (power-averaged) loudness
    loudnessISO3PowAvg = ((loudnessISO3TimeVar.iloc[int(start_skipT*sampleRateLoudISO3):
                                                    -int(end_skipT*sampleRateLoudISO3)]**(1/np.log10(2))).sum(axis=0)
                          / len(loudnessISO3TimeVar.iloc[int(start_skipT*sampleRateLoudISO3):
                                                          -int(end_skipT*sampleRateLoudISO3)]))**np.log10(2)
    # loudnessISO3PowAvg = loudnessISO3PowAvg.iloc[0]  # convert series to float

    # Calculate overall Aures+ISO532-3 sharpness from 2-channel time-varying
    # sharpness
    sharpAISO3TimeVar = pd.DataFrame(workbookdata['Sheet10'].iloc[13:, 1:3].values,
                                      columns=workbookdata['Sheet10'].iloc[12, 1:3],
                                      index=workbookdata['Sheet10'].iloc[13:, 0])
    # 2-channel overall (mean) sharpness
    sharpAISO3Mean = sharpAISO3TimeVar.iloc[int(start_skipT*sampleRateLoudISO3):
                                            -int(end_skipT*sampleRateLoudISO3)].mean(axis=0)
    # max of l/r channel overall (mean) sharpness
    sharpAISO3MeanMaxLR = sharpAISO3Mean.max()
    # 2-channel overall (power-averaged) sharpness
    sharpAISO3PowAvg = ((sharpAISO3TimeVar.iloc[int(start_skipT*sampleRateLoudISO3):
                                                -int(end_skipT*sampleRateLoudISO3)]**(1/np.log10(2))).sum(axis=0)
                        / len(sharpAISO3TimeVar.iloc[int(start_skipT*sampleRateLoudISO3):
                                                      -int(end_skipT*sampleRateLoudISO3)]))**np.log10(2)
    # max of l/r channel overall (power-averaged) sharpness
    sharpAISO3PowAvgMaxLR = sharpAISO3PowAvg.max()
    # 2-channel overall 5% exceeded sharpness
    sharpAISO305Ex = sharpAISO3TimeVar.iloc[int(start_skipT*sampleRateLoudISO3):
                                            -int(end_skipT*sampleRateLoudISO3)].quantile(q=0.95)
    # max of l/r channel overall 5% exceeded sharpness
    sharpAISO305ExMaxLR = sharpAISO305Ex.max()
    
    # Calculate overall Sottek Hearing Model impulsiveness from 2-channel
    # time-varying impulsiveness
    impulsiveHMSTimeVar = pd.DataFrame(workbookdata['Sheet12'].iloc[13:, 1:3].values,
                                        columns=workbookdata['Sheet12'].iloc[12, 1:3],
                                        index=workbookdata['Sheet12'].iloc[13:, 0])
    # 2-channel overall (power-averaged) impulsiveness
    impulsiveHMSPowAvg = ((impulsiveHMSTimeVar.iloc[int(np.ceil(sampleRateImpulsHMS*start_skipT)):
                                                    -int(np.ceil(sampleRateImpulsHMS*end_skipT))]**(1/np.log10(2))).sum(axis=0)
                          / len(impulsiveHMSTimeVar.iloc[int(np.ceil(sampleRateImpulsHMS*start_skipT)):
                                                          -int(np.ceil(sampleRateImpulsHMS*end_skipT))]))**np.log10(2)
    # max of l/r overall (power-averaged) impulsiveness
    impulsiveHMSPowAvgMaxLR = impulsiveHMSPowAvg.max()
    # 2-channel overall (mean) impulsiveness
    impulsiveHMSMean = impulsiveHMSTimeVar.iloc[int(np.ceil(sampleRateImpulsHMS*start_skipT)):
                                                -int(np.ceil(sampleRateImpulsHMS*end_skipT))].mean(axis=0)
    # max of l/r channel overall (mean) impulsiveness
    impulsiveHMSMeanMaxLR = impulsiveHMSMean.max()
    # max of l/r channel overall (max) impulsiveness
    impulsiveHMSMaxMaxLR = impulsiveHMSTimeVar.iloc[int(np.ceil(sampleRateImpulsHMS*start_skipT)):
                                                    -int(np.ceil(sampleRateImpulsHMS*end_skipT))].max(axis=None)
    # 2-channel overall 5% exceeded impulsiveness
    impulsiveHMS05Ex = impulsiveHMSTimeVar.iloc[int(np.ceil(sampleRateImpulsHMS*start_skipT)):
                                                -int(np.ceil(sampleRateImpulsHMS*end_skipT))].quantile(q=0.95)
    # max of l/r channel overall 5% exceeded impulsiveness
    impulsiveHMS05ExMaxLR = impulsiveHMS05Ex.max()

    # add results to output DataFrame
    dataByStim.loc[stimulus, 'LoudECMAHMSPowAvgBin'] = loudnessHMSPowAvg
    dataByStim.loc[stimulus, 'TonalECMAHMSAvgMaxLR'] = tonalityHMSAvgMaxLR
    dataByStim.loc[stimulus, 'TonalECMAHMS05ExMaxLR'] = tonalityHMS05ExMaxLR
    dataByStim.loc[stimulus, 'RoughECMAHMS10ExBin'] = roughHMS10Ex
    dataByStim.loc[stimulus, 'RoughECMAHMS05ExBin'] = roughHMS05Ex
    dataByStim.loc[stimulus, 'FluctstrHMS10ExBin'] = fluctStrHMSTime10Ex
    dataByStim.loc[stimulus, 'FluctstrHMS05ExBin'] = fluctStrHMSTime05Ex
    dataByStim.loc[stimulus, 'LoudISO105ExMaxLR'] = loudnessISO105ExMaxLR
    dataByStim.loc[stimulus, 'LoudISO1PowAvgMaxLR'] = loudnessISO1PowAvgMaxLR
    dataByStim.loc[stimulus, 'LoudISO3PowAvgMaxLR'] = loudnessISO3PowAvg
    dataByStim.loc[stimulus, 'SharpAuresISO3AvgMaxLR'] = sharpAISO3MeanMaxLR
    dataByStim.loc[stimulus, 'SharpAuresISO3PowAvgMaxLR'] = sharpAISO3PowAvgMaxLR
    dataByStim.loc[stimulus, 'SharpAuresISO305ExMaxLR'] = sharpAISO305ExMaxLR
    dataByStim.loc[stimulus, 'ImpulsHMSPowAvgMaxLR'] = impulsiveHMSPowAvgMaxLR
    dataByStim.loc[stimulus, 'ImpulsHMSAvgMaxLR'] = impulsiveHMSMeanMaxLR
    dataByStim.loc[stimulus, 'ImpulsHMSMaxMaxLR'] = impulsiveHMSMaxMaxLR
    dataByStim.loc[stimulus, 'ImpulsHMS05ExMaxLR'] = impulsiveHMS05ExMaxLR

    
    # calculation section for SQM differences
    # NOTE: THIS SECTION RELIES ON THE ALPHABETIC ORDER OF THE STIMULI FILES AS
    # ORIGINALLY NAMED: EACH AMBIENT SOUND FILE PRECEDING THE CORRESPONDING
    # COMBINED STIMULI FILES. IF THE ORDERING OR FILE NAMING IS CHANGED, THIS
    # CALCULATION WILL BE INVALIDATED AND *MAY ALSO* CAUSE AN ERROR.
    # a rolling 50 ms window is applied to average the SQM values over time -
    # this is to reduce uncertainty due to imperfect time-alignment between the
    # ambient vs combined stimuli files (all are from recordings, so there will
    # be some slippage due to imperfect editing)
    if stimulus in ["A1_CALBIN_Pa.wav", "A2_CALBIN_Pa.wav",
                    "B2_CALBIN_Pa.wav"]:
        # NOTE: we could dropna() the first <windowT values, but these will be
        # ignored anyway in the statistical analysis, assuming start_skipT >
        # windowT

        # calculate moving average values for ambient stimulus
        ambSpecTonalityHMSLMovAvg = specTonalityHMSL.T.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean().T
        ambSpecTonalityHMSRMovAvg = specTonalityHMSR.T.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean().T
        ambSpecRoughHMSBinMovAvg = specRoughHMSBin.T.rolling(window=int(np.ceil(sampleRateRoughECMA*windowT))).mean().T
        ambSpecFluctStrHMSLMovAvg = specFluctStrHMSL.T.rolling(window=int(np.ceil((sampleRateFluctHMS)*windowT))).mean().T
        ambSpecFluctStrHMSRMovAvg = specFluctStrHMSR.T.rolling(window=int(np.ceil((sampleRateFluctHMS)*windowT))).mean().T
        ambSharpAISO3TimeVarMovAvg = sharpAISO3TimeVar.rolling(window=int(np.ceil(sampleRateLoudISO3*windowT))).mean()
        ambImpulsiveHMSTimeVarMovAvg = impulsiveHMSTimeVar.rolling(window=int(np.ceil(sampleRateImpulsHMS*windowT))).mean()
        
    elif stimulus[0:3] in ["A1_", "A2_", "B2_"]:
        # calculation moving average values for combined stimulus
        specTonalityHMSLMovAvg = specTonalityHMSL.T.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean().T
        specTonalityHMSRMovAvg = specTonalityHMSR.T.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean().T
        specRoughHMSBinMovAvg = specRoughHMSBin.T.rolling(window=int(np.ceil(sampleRateRoughECMA*windowT))).mean().T
        specFluctStrHMSLMovAvg = specFluctStrHMSL.T.rolling(window=int(np.ceil((sampleRateFluctHMS)*windowT))).mean().T
        specFluctStrHMSRMovAvg = specFluctStrHMSR.T.rolling(window=int(np.ceil((sampleRateFluctHMS)*windowT))).mean().T
        sharpAISO3TimeVarMovAvg = sharpAISO3TimeVar.rolling(window=int(np.ceil(sampleRateLoudISO3*windowT))).mean()
        impulsiveHMSTimeVarMovAvg = impulsiveHMSTimeVar.rolling(window=int(np.ceil(sampleRateImpulsHMS*windowT))).mean()
        
        # calculate differences and make negative values 0
        dSpecTonalityHMSL = np.maximum(specTonalityHMSLMovAvg
                                        - ambSpecTonalityHMSLMovAvg, 0)
        dSpecTonalityHMSR = np.maximum(specTonalityHMSRMovAvg
                                        - ambSpecTonalityHMSRMovAvg, 0)
        dSpecRoughHMSBin = np.maximum(specRoughHMSBinMovAvg
                                      - ambSpecRoughHMSBinMovAvg, 0)
        dSpecFluctStrHMSL = np.maximum(specFluctStrHMSLMovAvg
                                        - ambSpecFluctStrHMSLMovAvg, 0)
        dSpecFluctStrHMSR = np.maximum(specFluctStrHMSRMovAvg
                                        - ambSpecFluctStrHMSRMovAvg, 0)
        dSharpAISO3TimeVar = np.maximum(sharpAISO3TimeVarMovAvg
                                        - ambSharpAISO3TimeVarMovAvg, 0)
        dImpulsiveHMSTimeVar = np.maximum(impulsiveHMSTimeVarMovAvg
                                          - ambImpulsiveHMSTimeVarMovAvg,
                                          0)

        # calculate aggregated difference values

        # 2-channel time-varing tonality
        dTonalityHMSTimeVar = pd.concat([dSpecTonalityHMSL.max(axis=0),
                                          dSpecTonalityHMSR.max(axis=0)],
                                        axis=1)
        # 2-channel time-averaged tonality (omitting T<=0.02)
        dTonalityHMSAvgL = dTonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                                    -int(np.ceil(sampleRateTonalECMA*end_skipT)),
                                                    0][dTonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                                                                -int(np.ceil(sampleRateTonalECMA*end_skipT)), 0]
                                                                                > 0.02].mean(axis=0)
        dTonalityHMSAvgR = dTonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                                    -int(np.ceil(sampleRateTonalECMA*end_skipT)),
                                                    1][dTonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                                                                -int(np.ceil(sampleRateTonalECMA*end_skipT)), 1]
                                                                                > 0.02].mean(axis=0)
        dTonalityHMSAvgMaxLR = max(dTonalityHMSAvgL, dTonalityHMSAvgR)
        # 2-channel 5% exceeded tonality (automatically omitting T<=0.02)
        dTonalityHMS05ExL = dTonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                                      -int(np.ceil(sampleRateTonalECMA*end_skipT)),
                                                      0][dTonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                                                                  -int(np.ceil(sampleRateTonalECMA*end_skipT)), 0]
                                                                                  > 0.02].quantile(q=0.95)
        dTonalityHMS05ExR = dTonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                                      -int(np.ceil(sampleRateTonalECMA*end_skipT)),
                                                      1][dTonalityHMSTimeVar.iloc[int(np.ceil(sampleRateTonalECMA*start_skipT)):
                                                                                  -int(np.ceil(sampleRateTonalECMA*end_skipT)), 1]
                                                                                  > 0.02].quantile(q=0.95)
        dTonalityHMS05ExMaxLR = max(dTonalityHMS05ExL, dTonalityHMS05ExR)

        # binaural time-varying roughness
        dRoughHMSTimeVar = dSpecRoughHMSBin.sum(axis=0)*0.5
        # binaural overall (90th percentile = 10% exceeded) roughness
        dRoughHMS10Ex = dRoughHMSTimeVar.iloc[int(np.ceil(sampleRateRoughECMA*start_skipT)):
                                              -int(np.ceil(sampleRateRoughECMA*end_skipT))].quantile(q=0.9)
        # binaural overall (95th percentile = 5% exceeded) roughness
        dRoughHMS05Ex = dRoughHMSTimeVar.iloc[int(np.ceil(sampleRateRoughECMA*start_skipT)):
                                              -int(np.ceil(sampleRateRoughECMA*end_skipT))].quantile(q=0.95)

        # binaural specific fluctuation strength
        # (using ECMA-418-2:2022 Equation 112 for roughness)
        dSpecFluctStrHMSBin = ((dSpecFluctStrHMSL**2
                                + dSpecFluctStrHMSR**2)/2).pow(1/2)
        # binaural time-dependent fluctuation strength
        # (using ECMA-418-2:2022 Equation 111 for roughness)
        dFluctStrHMSTimeVar = dSpecFluctStrHMSBin.sum(axis=0)*0.5
        # binaural overall (90th percentile = 10% exceeded) fluctuation strength
        # (using ECMA-418-2:2022 Section 7.1.8 for roughness)
        dFluctStrHMSTime10Ex = dFluctStrHMSTimeVar.iloc[int(np.ceil(start_skipT*(sampleRateFluctHMS))):
                                                        -int(np.ceil(end_skipT*(sampleRateFluctHMS)))].quantile(q=0.9)
        # binaural overall (95th percentile = 5% exceeded) fluctuation strength
        dFluctStrHMSTime05Ex = dFluctStrHMSTimeVar.iloc[int(np.ceil(start_skipT*(sampleRateFluctHMS))):
                                                        -int(np.ceil(end_skipT*(sampleRateFluctHMS)))].quantile(q=0.95)

        # 2-channel overall (mean) sharpness
        dSharpAISO3Mean = dSharpAISO3TimeVar.iloc[int(start_skipT*sampleRateLoudISO3):
                                                  -int(end_skipT*sampleRateLoudISO3)].mean(axis=0)
        # max of l/r channel overall (mean) sharpness
        dSharpAISO3MeanMaxLR = dSharpAISO3Mean.max()
        # 2-channel overall (power-averaged) sharpness
        dSharpAISO3PowAvg = ((dSharpAISO3TimeVar.iloc[int(start_skipT*sampleRateLoudISO3):
                                                      -int(end_skipT*sampleRateLoudISO3)]**(1/np.log10(2))).sum(axis=0)
                              / len(dSharpAISO3TimeVar.iloc[int(start_skipT*sampleRateLoudISO3):
                                                            -int(end_skipT*sampleRateLoudISO3)]))**np.log10(2)
        # max of l/r channel overall (power-averaged) sharpness
        dSharpAISO3PowAvgMaxLR = dSharpAISO3PowAvg.max()
        # 2-channel overall 5% exceeded sharpness
        dSharpAISO305Ex = dSharpAISO3TimeVar.iloc[int(start_skipT*sampleRateLoudISO3):
                                                  -int(end_skipT*sampleRateLoudISO3)].quantile(q=0.95)
        # max of l/r channel overall 5% exceeded sharpness
        dSharpAISO305ExMaxLR = dSharpAISO305Ex.max()

        # 2-channel overall (power-averaged) impulsiveness
        dImpulsiveHMSPowAvg = ((dImpulsiveHMSTimeVar.iloc[int(np.ceil(sampleRateImpulsHMS*start_skipT)):
                                                          -int(np.ceil(sampleRateImpulsHMS*end_skipT))]**(1/np.log10(2))).sum(axis=0)
                              / len(dImpulsiveHMSTimeVar.iloc[int(np.ceil(sampleRateImpulsHMS*start_skipT)):-int(np.ceil(sampleRateImpulsHMS*end_skipT))]))**np.log10(2)
        # max of l/r overall (power-averaged) impulsiveness
        dImpulsiveHMSPowAvgMaxLR = dImpulsiveHMSPowAvg.max()
        # 2-channel overall (mean) impulsiveness
        dImpulsiveHMSMean = dImpulsiveHMSTimeVar.iloc[int(np.ceil(sampleRateImpulsHMS*start_skipT)):
                                                      -int(np.ceil(sampleRateImpulsHMS*end_skipT))].mean(axis=0)
        # max of l/r channel overall (mean) impulsiveness
        dImpulsiveHMSMeanMaxLR = dImpulsiveHMSMean.max()
        # max of l/r channel overall (max) impulsiveness
        dImpulsiveHMSMaxMaxLR = dImpulsiveHMSTimeVar.iloc[int(np.ceil(sampleRateImpulsHMS*start_skipT)):
                                                          -int(np.ceil(sampleRateImpulsHMS*end_skipT))].max(axis=None)
        # 2-channel overall 5% exceeded impulsiveness
        dImpulsiveHMS05Ex = dImpulsiveHMSTimeVar.iloc[int(np.ceil(sampleRateImpulsHMS*start_skipT)):
                                                      -int(np.ceil(sampleRateImpulsHMS*end_skipT))].quantile(q=0.95)
        # max of l/r channel overall 5% exceeded impulsiveness
        dImpulsiveHMS05ExMaxLR = dImpulsiveHMS05Ex.max()

        # add results to output DataFrame
        dataByStim.loc[stimulus, 'dTonalECMAHMSAvgMaxLR'] = dTonalityHMSAvgMaxLR
        dataByStim.loc[stimulus, 'dTonalECMAHMS05ExMaxLR'] = dTonalityHMS05ExMaxLR
        dataByStim.loc[stimulus, 'dRoughECMAHMS10ExBin'] = dRoughHMS10Ex
        dataByStim.loc[stimulus, 'dRoughECMAHMS05ExBin'] = dRoughHMS05Ex
        dataByStim.loc[stimulus, 'dFluctstrHMS10ExBin'] = dFluctStrHMSTime10Ex
        dataByStim.loc[stimulus, 'dFluctstrHMS05ExBin'] = dFluctStrHMSTime05Ex
        dataByStim.loc[stimulus, 'dSharpAuresISO3AvgMaxLR'] = dSharpAISO3MeanMaxLR
        dataByStim.loc[stimulus, 'dSharpAuresISO3PowAvgMaxLR'] = dSharpAISO3PowAvgMaxLR
        dataByStim.loc[stimulus, 'dSharpAuresISO305ExMaxLR'] = dSharpAISO305ExMaxLR
        dataByStim.loc[stimulus, 'dImpulsHMSPowAvgMaxLR'] = dImpulsiveHMSPowAvgMaxLR
        dataByStim.loc[stimulus, 'dImpulsHMSAvgMaxLR'] = dImpulsiveHMSMeanMaxLR
        dataByStim.loc[stimulus, 'dImpulsHMSMaxMaxLR'] = dImpulsiveHMSMaxMaxLR
        dataByStim.loc[stimulus, 'dImpulsHMS05ExMaxLR'] = dImpulsiveHMS05ExMaxLR

# end of for loop over psychoacoustic metrics xlsx files

# add UAS-only and ambient-only data to combined stimuli
# ------------------------------------------------------
indicesDiffPsycho = ["dTonalECMAHMSAvgMaxLR",
                     "dTonalECMAHMS05ExMaxLR",
                     "dRoughECMAHMS10ExBin",
                     "dRoughECMAHMS05ExBin",
                     "dFluctstrHMS10ExBin",
                     "dFluctstrHMS05ExBin",
                     "dSharpAuresISO3AvgMaxLR",
                     "dSharpAuresISO3PowAvgMaxLR",
                     "dSharpAuresISO305ExMaxLR",
                     "dImpulsHMSPowAvgMaxLR",
                     "dImpulsHMSAvgMaxLR",
                     "dImpulsHMSMaxMaxLR",
                     "dImpulsHMS05ExMaxLR"]

indicesAbsPsycho = [index for index in indicesPsycho
                    if index not in indicesDiffPsycho]

indices = indicesAcoustic + indicesAbsPsycho

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
UASAmb = pd.concat([UASonly, Ambonly], axis=1)

# insert zeros for No UAS stimuli UAS SQMs
NoUASStims = ['A1_CALBIN_Pa.wav', 'A2_CALBIN_Pa.wav', 'B2_CALBIN_Pa.wav']
indicesUASPsycho = ["UAS" + index for index in indicesPsycho if index not in indicesDiffPsycho]
UASAmb.loc[NoUASStims, indicesUASPsycho] = 0

# merge into output
dataByStim = dataByStim.merge(UASAmb.astype(float), how='outer',
                              left_index=True, right_index=True)

# insert zeros for SQM differences
dataByStim.loc[NoUASStims, indicesDiffPsycho] = 0

# calculate level differences and ratios between UAS and ambient
dataByStim['UASLAeqdiffAmbLAF90'] = dataByStim['UASLAeqMaxLR'] - dataByStim['AmbLAF90ExMaxLR']
dataByStim['UASLASmaxdiffAmbLAF90'] = dataByStim['UASLASmaxMaxLR'] - dataByStim['AmbLAF90ExMaxLR']
dataByStim['UASLASmaxdiffAmbLAF50'] = dataByStim['UASLASmaxMaxLR'] - dataByStim['AmbLAF50ExMaxLR']
dataByStim['UASLASmaxdiffAmbLAeq'] = dataByStim['UASLASmaxMaxLR'] - dataByStim['AmbLAeqMaxLR']
dataByStim['EPNLdiff'] = dataByStim['UASEPNLMaxLR'] - dataByStim['AmbEPNLMaxLR']
dataByStim['EPNLNoTdiff'] = dataByStim['UASEPNLNoTMaxLR'] - dataByStim['AmbEPNLNoTMaxLR']
dataByStim['PNLTmaxdiff'] = dataByStim['UASPNLTmaxMaxLR'] - dataByStim['AmbPNLTmaxMaxLR']


# ---------------------
# Partial loudness data
# ---------------------

# open csv file selection dialog and assign filepaths to list
# PROJECT NOTE: the calculated files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Analysis/deeuu_loudness/output?csf=1&web=1&e=ZvblMt
app = QApplication(sys.argv)
fileExts = "*TermPartialLoudness.pkl"
filelist = list(QFileDialog.getOpenFileNames(filter=fileExts))[0]
filelist.sort()

filenames = [filepath.split('/')[-1] for filepath in filelist]
filenamesST = [filename for filename in filenames
               if filename.find("Short") != -1]
filenamesLT = [filename for filename in filenames
               if filename.find("Long") != -1]

# setup results DataFrame
partLoudnessST = pd.DataFrame(index=filenamesST,
                              columns=['UASPartLoudGMSTPowAvg',
                                       'UASPartLoudGMST05Ex'])
partLoudnessLT = pd.DataFrame(index=filenamesLT,
                              columns=['UASPartLoudGMLTPowAvg',
                                       'UASPartLoudGMLT05Ex'])
# the deeuu Glasberg&Moore loudness outputs are all sampled at 1 ms
# in each of the deeuu loudness outputs, the first 0.031 seconds is
# redundant 'negative time' and the final 0.169 seconds is redundant added time
# skip indices
skip_starti = 31
skip_endi = -169

sampleRatePartLoudGM = 1e3

# loop over files to extract data
for ii, file in enumerate(filelist):
    
    if filenames[ii].find("A") == 0:
        timeTotal = 25.0
    elif filenames[ii].find("B") == 0:
        timeTotal = 75.0

    partN = pd.read_pickle(file)
    partN = pd.Series(data=partN[skip_starti:skip_endi, 1],
                      index=np.arange(0, timeTotal, 1e-3))

    partialLoudnessPowAvg = (np.sum(partN.iloc[int(start_skipT*sampleRatePartLoudGM)
                                               :int(-end_skipT*sampleRatePartLoudGM)]**(1/np.log10(2)),
                                    axis=0)
                             / len(partN.iloc[int(start_skipT*sampleRatePartLoudGM)
                                              :int(-end_skipT*sampleRatePartLoudGM)]))**np.log10(2)
    
    partialLoudness05Ex = partN.iloc[int(start_skipT*sampleRatePartLoudGM)
                                     :int(-end_skipT*sampleRatePartLoudGM)].quantile(q=0.95)

    if file.find("Short") != -1:
        partLoudnessST.loc[filenames[ii], 'UASPartLoudGMSTPowAvg'] = partialLoudnessPowAvg
        partLoudnessST.loc[filenames[ii], 'UASPartLoudGMST05Ex'] = partialLoudness05Ex
    elif file.find("Long") != -1:
        partLoudnessLT.loc[filenames[ii], 'UASPartLoudGMLTPowAvg'] = partialLoudnessPowAvg
        partLoudnessLT.loc[filenames[ii], 'UASPartLoudGMLT05Ex'] = partialLoudness05Ex

# reindex partial loudness DataFrame to match recording files and merge into
# output DataFrame
partLoudnessST['newindex'] = [file.replace("ShortTermPartialLoudness.pkl", ".wav")
                              for file in list(partLoudnessST.index)]
partLoudnessST.set_index('newindex', inplace=True)
partLoudnessLT['newindex'] = [file.replace("LongTermPartialLoudness.pkl", ".wav")
                              for file in list(partLoudnessLT.index)]
partLoudnessLT.set_index('newindex', inplace=True)

partLoudness = partLoudnessST.merge(right=partLoudnessLT, how='outer',
                                    left_index=True, right_index=True)

dataByStim = dataByStim.merge(partLoudness.astype(float), how='outer',
                              left_index=True, right_index=True)

# insert zeros for No UAS stimuli partial loudness
dataByStim.loc[NoUASStims, ['UASPartLoudGMSTPowAvg', 'UASPartLoudGMST05Ex',
                            'UASPartLoudGMLTPowAvg', 'UASPartLoudGMLT05Ex'] ] = 0

# move SQM difference metrics to after level difference metrics
indicesDiffPsycho.reverse()
for col in indicesDiffPsycho:
    popcol = dataByStim.pop(col)
    dataByStim.insert(loc=(dataByStim.columns.get_loc(partLoudness.columns[-1])
                           + 1),
                      column=col, value=popcol)
indicesDiffPsycho.reverse()  # revert list for use later   


# --------------------------------
# Detectability index calculations
# --------------------------------

detectEfficiencyData = pd.DataFrame(data=np.array([[32.3, 39.8, 51.2, 66.2,
                                                   89.1, 119.9, 161.4, 217.2,
                                                   292.3, 393.4, 529.5, 654.1,
                                                   800.5, 979.9, 1290.8,
                                                   1737.3, 2338.2, 3146.9,
                                                   4235.4, 5700.3, 7671.9,
                                                   10325.5, 13441.1],
                                                  [0.134, 0.163, 0.203, 0.247,
                                                   0.293, 0.334, 0.369, 0.394,
                                                   0.411, 0.422, 0.429, 0.430,
                                                   0.431, 0.430, 0.424, 0.414,
                                                   0.397, 0.373, 0.343, 0.310,
                                                   0.276, 0.243, 0.213]]).T,
                                    columns=["Hz", "eta"])

def detectEfficiencyFunc(freq, a, b, c, d, e, f):
    return (freq/800)**a/(b + c*(freq/800)**d)**e + f

p0 = [19, 0.02, 1, 0.14, 140, -0.9]

fm, f1, f2 = dsp.noct.noctf(20, 20e3, 3)

detectEffCurveFit, _ = optimize.curve_fit(f=detectEfficiencyFunc,
                                          xdata=detectEfficiencyData['Hz'],
                                          ydata=detectEfficiencyData['eta'],
                                          p0=p0, maxfev=1000000,
                                          bounds=(-10, 10000))
detectEffCurveTestFit = detectEfficiencyFunc(detectEfficiencyData['Hz'],
                                             detectEffCurveFit[0],
                                             detectEffCurveFit[1],
                                             detectEffCurveFit[2],
                                             detectEffCurveFit[3],
                                             detectEffCurveFit[4],
                                             detectEffCurveFit[5])



# -------------
# Response data
# -------------

# open xlsx file selection dialog and assign filepath
# PROJECT NOTE: the results files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Test_files/Results?csf=1&web=1&e=f3ejcD
app = QApplication(sys.argv)
fileExts = "*.xlsx"
filepath = QFileDialog.getOpenFileName(filter=fileExts,
                                       caption=r"Select test response data files")[0]


# Part A
# ------

# read in data and add column indicating the stimulus recording file
partAResponses = pd.read_excel(io=filepath, sheet_name="PtAResponse",
                               header=0)
partAResponses['Recording'] = [StimFile.split('.')[0] + "_CALBIN_Pa.wav"
                               for StimFile in partAResponses['Stim File']]

partAData = partAResponses.drop(columns=['Typology', 'Aircraft', 'Other'])
partAData.columns = [s.replace(' ', '') for s in partAData.columns]

# add highly annoyed data
partAData['HighAnnoy'] = (partAData['Annoyance'] >= 8).astype(int)
partAData.insert(loc=partAData.columns.get_loc('UAS_noticed'),
                 column='HighAnnoy', value=partAData.pop('HighAnnoy'))

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
    partAHighAnnoy = (partAAnnoy >= 8).astype(int)
    partAHighAnnoy.columns = partAHighAnnoy.columns.str.replace("Annoyance",
                                                                "HighAnnoy")
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
    highAnnoyAgg = pd.DataFrame(data=[[np.sum(partAHighAnnoy.values,
                                              axis=1)[0],
                                       np.mean(partAHighAnnoy.values,
                                               axis=1)[0]]],
                                columns=['HighAnnoyTotal',
                                         'HighAnnoyProportion'],
                                index=[file])
    noticeAgg = pd.DataFrame(data=[[np.sum(partANotice.values, axis=1)[0],
                                    np.mean(partANotice.values, axis=1)[0]]],
                             columns=['NoticedTotal', 'NoticedProportion'],
                             index=[file])

    # add results to DataFrame
    if ii == 0:
        partA = partA.join([partAValence, partAArousal, partAAnnoy,
                            partAHighAnnoy,
                            partANotice, valenceAgg, arousalAgg,
                            annoyAgg, highAnnoyAgg, noticeAgg])
        partAstats = partAstats.join([valenceSWtest, arousalSWtest,
                                      annoySWtest, valenceAgg, arousalAgg,
                                      annoyAgg, highAnnoyAgg, noticeAgg])
    else:
        partAstats.loc[file, valenceSWtest.columns] = valenceSWtest.loc[file]
        partAstats.loc[file, arousalSWtest.columns] = arousalSWtest.loc[file]
        partAstats.loc[file, annoySWtest.columns] = annoySWtest.loc[file]
        partAstats.loc[file, valenceAgg.columns] = valenceAgg.loc[file]
        partAstats.loc[file, arousalAgg.columns] = arousalAgg.loc[file]
        partAstats.loc[file, annoyAgg.columns] = annoyAgg.loc[file]
        partAstats.loc[file, highAnnoyAgg.columns] = highAnnoyAgg.loc[file]
        partAstats.loc[file, noticeAgg.columns] = noticeAgg.loc[file]

        partA.loc[file, partAValence.columns] = partAValence.loc[file]
        partA.loc[file, partAArousal.columns] = partAArousal.loc[file]
        partA.loc[file, partAAnnoy.columns] = partAAnnoy.loc[file]
        partA.loc[file, partAHighAnnoy.columns] = partAHighAnnoy.loc[file]
        partA.loc[file, partANotice.columns] = partANotice.loc[file]
        partA.loc[file, valenceAgg.columns] = valenceAgg.loc[file]
        partA.loc[file, arousalAgg.columns] = arousalAgg.loc[file]
        partA.loc[file, annoyAgg.columns] = annoyAgg.loc[file]
        partA.loc[file, highAnnoyAgg.columns] = highAnnoyAgg.loc[file]
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

# add highly annoyed data
partBData['HighAnnoy'] = (partBData['Annoyance'] >= 8).astype(int)
partBData.insert(loc=partBData.columns.get_loc('Recording'),
                 column='HighAnnoy', value=partBData.pop('HighAnnoy'))

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
    partBHighAnnoy = (partBAnnoy >= 8).astype(int)
    partBHighAnnoy.columns = partBHighAnnoy.columns.str.replace("Annoyance",
                                                                "HighAnnoy")
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
    highAnnoyAgg = pd.DataFrame(data=[[np.sum(partBHighAnnoy.values,
                                              axis=1)[0],
                                       np.mean(partBHighAnnoy.values,
                                               axis=1)[0]]],
                                columns=['HighAnnoyTotal',
                                         'HighAnnoyProportion'],
                                index=[file])

    # add results to DataFrame
    if ii == 0:
        partB = partB.join([partBValence, partBArousal, partBAnnoy,
                            partBHighAnnoy,
                            valenceAgg, arousalAgg, annoyAgg, highAnnoyAgg])
        partBstats = partBstats.join([valenceSWtest, arousalSWtest,
                                      annoySWtest, valenceAgg, arousalAgg,
                                      annoyAgg, highAnnoyAgg])
    else:
        partBstats.loc[file, valenceSWtest.columns] = valenceSWtest.loc[file]
        partBstats.loc[file, arousalSWtest.columns] = arousalSWtest.loc[file]
        partBstats.loc[file, annoySWtest.columns] = annoySWtest.loc[file]
        partBstats.loc[file, valenceAgg.columns] = valenceAgg.loc[file]
        partBstats.loc[file, arousalAgg.columns] = arousalAgg.loc[file]
        partBstats.loc[file, annoyAgg.columns] = annoyAgg.loc[file]
        partBstats.loc[file, highAnnoyAgg.columns] = highAnnoyAgg.loc[file]
        
        partB.loc[file, partBValence.columns] = partBValence.loc[file]
        partB.loc[file, partBArousal.columns] = partBArousal.loc[file]
        partB.loc[file, partBAnnoy.columns] = partBAnnoy.loc[file]
        partB.loc[file, partBHighAnnoy.columns] = partBHighAnnoy.loc[file]
        partB.loc[file, valenceAgg.columns] = valenceAgg.loc[file]
        partB.loc[file, arousalAgg.columns] = arousalAgg.loc[file]
        partB.loc[file, annoyAgg.columns] = annoyAgg.loc[file]
        partB.loc[file, highAnnoyAgg.columns] = highAnnoyAgg.loc[file]
    
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
# replace problem columns names
postTestResponses.columns = postTestResponses.columns.str.replace(" ", "_")

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
dataByStim.loc[dataByStim.index.str.find("_H520_") != -1, 'UASType'] = "H520"
dataByStim.loc[dataByStim.index.str.find("_M300_") != -1, 'UASType'] = "M300"
dataByStim.loc[dataByStim.index.str.find("_T150_") != -1, 'UASType'] = "T150"
dataByStim.loc[np.logical_or(dataByStim.index.str.find("B2_CALBIN") != -1,
                             np.logical_or(dataByStim.index.str.find("A1_CALBIN") != -1,
                                           dataByStim.index.str.find("A2_CALBIN") != -1)),
               'UASType'] = "No UAS"
dataByStim['UASType'] = pd.Categorical(dataByStim['UASType'], ["No UAS",
                                                               "H520",
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


# merge response and stimuli data into 'by participant' test datasets, and
# save to file

partADataBySubj = pd.merge(left=partAData,
                           right=dataByStimTestA.loc[:,
                                                     :indicesDiffPsycho[-1]],
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
                                                     :indicesDiffPsycho[-1]],
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

# Re-save filtered datasets for noticeability analysis only
omitParticipants = [2, 5, 21, 32, 34, 36, 38, 39, 40, 44, 45]
omitColumns = ["UAS_noticed_" + str(partID) for partID in omitParticipants]

partADataBySubjNotice = partADataBySubj.loc[~partADataBySubj['ID#'].isin(omitParticipants), :]
partADataBySubjNoticeFilePath = os.path.join(outFilePath,
                                             "refmap_listest1_testdataANoticeFilt_BySubj.csv")
partADataBySubjNotice.to_csv(partADataBySubjNoticeFilePath, index=False)

omitColumns = omitColumns + (["Arousal_" + str(partID) for partID in omitParticipants])
omitColumns = omitColumns + (["Valence_" + str(partID) for partID in omitParticipants])
omitColumns = omitColumns + (["Annoyance_" + str(partID) for partID in omitParticipants])
omitColumns = omitColumns + (["HighAnnoy_" + str(partID) for partID in omitParticipants])
dataByStimTestANotice = dataByStimTestA.drop(labels=omitColumns, axis=1)
dataByStimTestANotice.drop(labels=['ArousalMean', 'ArousalMedian',
                                   'ValenceMean', 'ValenceMedian',
                                   'AnnoyMean', 'AnnoyMedian',
                                   'HighAnnoyTotal', 'HighAnnoyProportion',
                                   'NoticedTotal', 'NoticedProportion'], axis=1,
                           inplace=True)

keepColumns = [col for col in dataByStimTestA.columns[dataByStimTestA.columns.str.find("UAS_noticed_") == 0]
               if col not in omitColumns]
dataByStimTestANotice['NoticedTotalFilt'] = dataByStimTestANotice[keepColumns].sum(axis=1)
dataByStimTestANotice['NoticedPropFilt'] = dataByStimTestANotice[keepColumns].mean(axis=1)
keepParticipants = [label.replace("UAS_noticed_", "") for label in keepColumns]
keepColumns = [label.replace("UAS_noticed_", "Arousal_") for label in keepColumns]
dataByStimTestANotice['ArousalMeanFilt'] = dataByStimTestANotice[keepColumns].mean(axis=1)
dataByStimTestANotice['ArousalMedianFilt'] = dataByStimTestANotice[keepColumns].median(axis=1)
keepColumns = [label.replace("Arousal_", "Valence_") for label in keepColumns]
dataByStimTestANotice['ValenceMeanFilt'] = dataByStimTestANotice[keepColumns].mean(axis=1)
dataByStimTestANotice['ValenceMedianFilt'] = dataByStimTestANotice[keepColumns].median(axis=1)
keepColumns = [label.replace("Valence_", "Annoyance_") for label in keepColumns]
dataByStimTestANotice['AnnoyMeanFilt'] = dataByStimTestANotice[keepColumns].mean(axis=1)
dataByStimTestANotice['AnnoyMedianFilt'] = dataByStimTestANotice[keepColumns].median(axis=1)
keepColumns = [label.replace("Annoyance_", "HighAnnoy_") for label in keepColumns]
dataByStimTestANotice['HighAnnoyTotalFilt'] = dataByStimTestANotice[keepColumns].sum(axis=1)
dataByStimTestANotice['HighAnnoyPropFilt'] = dataByStimTestANotice['HighAnnoyTotalFilt']/len(keepParticipants)

dataByStimTestANoticeFilePath = os.path.join(outFilePath,
                                             "refmap_listest1_testdataANoticeFilt_ByStim.csv")
dataByStimTestANotice.to_csv(dataByStimTestANoticeFilePath)

# save pre-test and post-test response data to separate files

preTestDataFilePath = os.path.join(outFilePath,
                                   "refmap_listest1_pretestdata.csv")
preTestResponses.to_csv(preTestDataFilePath, index=False)

postTestDataFilePath = os.path.join(outFilePath,
                                   "refmap_listest1_posttestdata.csv")
postTestResponses.to_csv(postTestDataFilePath, index=False)

# form wide format datasets for each part and outcome
# Part A
partAAnnoyDataBySubjWide = partADataBySubj.pivot(index='ID#',
                                                 columns='Recording',
                                                 values='Annoyance')
partAValenceDataBySubjWide = partADataBySubj.pivot(index='ID#',
                                                   columns='Recording',
                                                   values='Valence')
partAArousalDataBySubjWide = partADataBySubj.pivot(index='ID#',
                                                   columns='Recording',
                                                   values='Arousal')
partANoticeDataBySubjWide = partADataBySubj.pivot(index='ID#',
                                                  columns='Recording',
                                                  values='UAS_noticed')
partANoticeFiltDataBySubjWide = partADataBySubjNotice.pivot(index='ID#',
                                                        columns='Recording',
                                                        values='UAS_noticed')

# merge all together for multivariate analysis
partADataBySubjWide = partADataBySubj.pivot(index='ID#',
                                            columns='Recording',
                                            values='Annoyance')
partADataBySubjWide.columns = partADataBySubjWide.columns.str.replace("_CALBIN_Pa.wav", "_annoy")
partADataBySubjWide = partADataBySubjWide.merge(partAValenceDataBySubjWide,
                                                on='ID#')
partADataBySubjWide.columns = partADataBySubjWide.columns.str.replace("_CALBIN_Pa.wav", "_valence")
partADataBySubjWide = partADataBySubjWide.merge(partAArousalDataBySubjWide,
                                                on='ID#')
partADataBySubjWide.columns = partADataBySubjWide.columns.str.replace("_CALBIN_Pa.wav", "_arousal")
partADataBySubjWide = partADataBySubjWide.merge(partANoticeDataBySubjWide,
                                                on='ID#')
partADataBySubjWide.columns = partADataBySubjWide.columns.str.replace("_CALBIN_Pa.wav", "_notice")

# merge with participant info and then save
partAAnnoyDataBySubjWide = partAAnnoyDataBySubjWide.merge(prePostTestResponses,
                                                          on='ID#')
partAValenceDataBySubjWide = partAValenceDataBySubjWide.merge(prePostTestResponses,
                                                              on='ID#')
partAArousalDataBySubjWide = partAArousalDataBySubjWide.merge(prePostTestResponses,
                                                              on='ID#')
partANoticeDataBySubjWide = partANoticeDataBySubjWide.merge(prePostTestResponses,
                                                            on='ID#')
partANoticeFiltDataBySubjWide = partANoticeFiltDataBySubjWide.merge(prePostTestResponses,
                                                                    on='ID#')
partADataBySubjWide = partADataBySubjWide.merge(prePostTestResponses,
                                                on='ID#')

partAAnnoyDataBySubjWide.to_csv(os.path.join(outFilePath,
                                             "refmap_listest1_AnnoyDataABySubjWide.csv"),
                                             index=False)
partAValenceDataBySubjWide.to_csv(os.path.join(outFilePath,
                                               "refmap_listest1_ValenceDataABySubjWide.csv"),
                                  index=False)
partAArousalDataBySubjWide.to_csv(os.path.join(outFilePath,
                                               "refmap_listest1_ArousalDataABySubjWide.csv"),
                                  index=False)
partANoticeDataBySubjWide.to_csv(os.path.join(outFilePath,
                                              "refmap_listest1_NoticeDataABySubjWide.csv"),
                                 index=False)
partANoticeFiltDataBySubjWide.to_csv(os.path.join(outFilePath,
                                                  "refmap_listest1_NoticeFiltDataABySubjWide.csv"),
                                     index=False)
partADataBySubjWide.to_csv(os.path.join(outFilePath,
                                        "refmap_listest1_DataABySubjWide.csv"),
                           index=False)

# filter for anomalous notice
partADataBySubjSubsetFiltWide = partADataBySubjWide.loc[~partADataBySubjWide.index.isin(omitParticipants), :]
partADataBySubjWide.to_csv(os.path.join(outFilePath,
                                        "refmap_listest1_DataABySubjSubsetFiltWide.csv"),
                           index=False)


# Part B
partBAnnoyDataBySubjWide = partBDataBySubj.pivot(index='ID#',
                                                 columns='Recording',
                                                 values='Annoyance')
partBValenceDataBySubjWide = partBDataBySubj.pivot(index='ID#',
                                                   columns='Recording',
                                                   values='Valence')
partBArousalDataBySubjWide = partBDataBySubj.pivot(index='ID#',
                                                   columns='Recording',
                                                   values='Arousal')

# merge all together for multivariate analysis
partBDataBySubjWide = partBDataBySubj.pivot(index='ID#',
                                            columns='Recording',
                                            values='Annoyance')
partBDataBySubjWide.columns = partBDataBySubjWide.columns.str.replace("_CALBIN_Pa.wav", "_annoy")
partBDataBySubjWide = partBDataBySubjWide.merge(partBValenceDataBySubjWide,
                                                on='ID#')
partBDataBySubjWide.columns = partBDataBySubjWide.columns.str.replace("_CALBIN_Pa.wav", "_valence")
partBDataBySubjWide = partBDataBySubjWide.merge(partBArousalDataBySubjWide,
                                                on='ID#')
partBDataBySubjWide.columns = partBDataBySubjWide.columns.str.replace("_CALBIN_Pa.wav", "_arousal")

# merge with participant info and then save
partBAnnoyDataBySubjWide = partBAnnoyDataBySubjWide.merge(prePostTestResponses,
                                                          on='ID#')
partBValenceDataBySubjWide = partBValenceDataBySubjWide.merge(prePostTestResponses,
                                                              on='ID#')
partBArousalDataBySubjWide = partBArousalDataBySubjWide.merge(prePostTestResponses,
                                                              on='ID#')
partBDataBySubjWide = partBDataBySubjWide.merge(prePostTestResponses,
                                                on='ID#')

partBAnnoyDataBySubjWide.to_csv(os.path.join(outFilePath,
                                             "refmap_listest1_AnnoyDataBBySubjWide.csv"),
                                index=False)
partBValenceDataBySubjWide.to_csv(os.path.join(outFilePath,
                                               "refmap_listest1_ValenceDataBBySubjWide.csv"),
                                  index=False)
partBArousalDataBySubjWide.to_csv(os.path.join(outFilePath,
                                               "refmap_listest1_ArousalDataBBySubjWide.csv"),
                                  index=False)
partBDataBySubjWide.to_csv(os.path.join(outFilePath,
                                        "refmap_listest1_DataBBySubjWide.csv"),
                           index=False)