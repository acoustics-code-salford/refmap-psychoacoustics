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
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

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

    # skip samples
    start_skips = int(start_skipT*sampleRatein)
    end_skips = int(end_skipT*sampleRatein)

    # get time vector
    dT = 1/sampleRatein
    endT = len(signal)/sampleRatein
    timeVector = np.arange(0, endT, dT)
    timeVectorSkip = timeVector[start_skips:-end_skips]

    # apply weighting filters
    signalA = dsp.filterFuncs.A_weight_T(signal, sampleRatein)
    signalmagAF = dsp.filterFuncs.time_weight(signalA, sampleRatein, tau=0.125)
    signalmagAS = dsp.filterFuncs.time_weight(signalA, sampleRatein, tau=1)

    # calculate weighted dB time series
    signaldBAF = 20*np.log10(signalmagAF[start_skips:-end_skips]/2e-5)
    signaldBAS = 20*np.log10(signalmagAS[start_skips:-end_skips]/2e-5)

    # calculate percentile metrics
    signalLAFmax = signaldBAF.max(axis=0)
    signalLAF5 = np.percentile(signaldBAF, q=95, axis=0)
    signalLAF10 = np.percentile(signaldBAF, q=90, axis=0)
    signalLAF25 = np.percentile(signaldBAF, q=75, axis=0)
    signalLAF50 = np.percentile(signaldBAF, q=50, axis=0)
    signalLAF75 = np.percentile(signaldBAF, q=25, axis=0)
    signalLAF90 = np.percentile(signaldBAF, q=10, axis=0)
    signalLAF95 = np.percentile(signaldBAF, q=5, axis=0)
    signalLASmax = signaldBAS.max(axis=0)

    # calculate energy metrics
    signalLAeq = 20*np.log10(np.sqrt((signalA[start_skips:
                                              -end_skips]**2).mean(axis=0))/2e-5)
    signalLAE = (signalLAeq
                 + 10*np.log10(len(signalA[start_skips:
                                           -end_skips])/sampleRatein))

    # ICAO LAE metrics
    signalLAFmaxLoc = signaldBAF.argmax(axis=0)
    signalLASmaxLoc = signaldBAS.argmax(axis=0)
    SlowLAEThresh = signalLASmax - 10
    FastLAEThresh = signalLAFmax - 10
    
    SlowIntStartL = timeVectorSkip[signaldBAS[:, 0] >= SlowLAEThresh[0]].min()
    SlowIntEndL = timeVectorSkip[signaldBAS[:, 0] >= SlowLAEThresh[0]].max()
    SlowIntStartR = timeVectorSkip[signaldBAS[:, 1] >= SlowLAEThresh[1]].min()
    SlowIntEndR = timeVectorSkip[signaldBAS[:, 1] >= SlowLAEThresh[1]].max()
    
    SlowIntTime = np.array([SlowIntEndL - SlowIntStartL,
                            SlowIntEndR - SlowIntStartR])
    
    # note that the indices here refer to the full signal without truncation (start_skips is not omitted from the sample index)
    SlowIntIndices = np.array([[int(SlowIntStartL*sampleRatein),
                                int(SlowIntStartR*sampleRatein)],
                               [int(SlowIntEndL*sampleRatein), 
                                int(SlowIntEndR*sampleRatein)]])
    
    signalICAOLAES = np.array([20*np.log10(np.sqrt(SlowIntTime[0]*(signalA[SlowIntIndices[0, 0]:
                                                                           SlowIntIndices[1, 0], 0]**2).mean(axis=0))/2e-5),
                                20*np.log10(np.sqrt(SlowIntTime[1]*(signalA[SlowIntIndices[0, 1]:
                                                                            SlowIntIndices[1, 1], 1]**2).mean(axis=0))/2e-5)])
    
    FastIntStartL = timeVectorSkip[signaldBAF[:, 0] >= FastLAEThresh[0]].min()
    FastIntEndL = timeVectorSkip[signaldBAF[:, 0] >= FastLAEThresh[0]].max()
    FastIntStartR = timeVectorSkip[signaldBAF[:, 1] >= FastLAEThresh[1]].min()
    FastIntEndR = timeVectorSkip[signaldBAF[:, 1] >= FastLAEThresh[1]].max()
    
    FastIntTime = np.array([FastIntEndL - FastIntStartL,
                            FastIntEndR - FastIntStartR])
    
    # note that the indices here refer to the full signal without truncation (start_skips is not omitted from the sample index)
    FastIntIndices = np.array([[int(FastIntStartL*sampleRatein),
                                int(FastIntStartR*sampleRatein)],
                               [int(FastIntEndL*sampleRatein), 
                                int(FastIntEndR*sampleRatein)]])
    
    signalTradLAEF = np.array([20*np.log10(np.sqrt(FastIntTime[0]*(signalA[FastIntIndices[0, 0]:
                                                                           FastIntIndices[1, 0], 0]**2).mean(axis=0))/2e-5),
                                20*np.log10(np.sqrt(FastIntTime[1]*(signalA[FastIntIndices[0, 1]:
                                                                            FastIntIndices[1, 1], 1]**2).mean(axis=0))/2e-5)])

    dataByStim.loc[filenames[ii], 'LAESICAOMaxLR'] = signalICAOLAES.max()
    dataByStim.loc[filenames[ii], 'LAEFTradMaxLR'] = signalTradLAEF.max()

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


    # calculate intermittency ratio for test sounds
    if (filenames[ii].find("A1") != -1 or filenames[ii].find("A2") != -1
        or filenames[ii].find("B2") != -1):
        # 1-second LAeq
        signalLAeq1s = 20*np.log10(np.sqrt(pd.DataFrame((signalA[start_skips:-end_skips]**2)).rolling(window=sampleRatein, step=sampleRatein).mean())/2e-5)
        intermitRatioK = signalLAeq + 3  # IR threshold
        mask = signalLAeq1s > intermitRatioK
        # calculate average intensity for events exceeding threshold
        IAeqEvents = np.nansum(10**(signalLAeq1s[mask]/10), axis=0)/(len(signalA[start_skips:-end_skips])/sampleRatein)
        # convert to sound pressure while avoiding log(0)
        LAeqEvents = 10*np.log10(IAeqEvents, out=np.zeros_like(IAeqEvents, dtype=np.float64),
                                 where=(IAeqEvents!=0))
        intermitRatio = 100*10**((LAeqEvents - signalLAeq)/10)
        intermitRatio[IAeqEvents == 0] = 0  # replace tiny decimal with 0
        dataByStim.loc[filenames[ii], 'IntermitRatioMaxLR'] = intermitRatio.max()
    else:
        dataByStim.loc[filenames[ii], 'IntermitRatioMaxLR'] = np.nan


# end of for loop over signal wav files


# ------------------------------------------------
# Psychoacoustic metrics single values calculation
# ------------------------------------------------

# output variables
indicesPsycho = ["LoudECMAHMSPowAvgBin",
                 "LoudISO105ExMaxLR",
                 "LoudISO1PowAvgMaxLR",
                 "LoudISO3PowAvgBin",
                 "TonalECMAHMSAvgMaxLR",
                 "TonalECMAHMS05ExMaxLR",
                 "TonalIntAvgMaxLR",
                 "TonalInt05ExMaxLR",
                 "TonLdECMAHMSPowAvgBin",
                 "TonLdECMAHMS05ExBin",
                 "RoughECMAHMS10ExBin",
                 "RoughECMAHMS05ExBin",
                 "FluctHMS10ExBin",
                 "FluctHMS05ExBin",
                 "SharpAuresISO3PowAvgMaxLR",
                 "SharpAuresISO305ExMaxLR",
                 "ImpulsHMSPowAvgMaxLR",
                 "ImpulsHMSAvgMaxLR",
                 "ImpulsHMS05ExMaxLR",
                 "dTonalECMAHMSAvgMaxLR",
                 "dTonalECMAHMS05ExMaxLR",
                 "dTonalIntAvgMaxLR",
                 "dTonalInt05ExMaxLR",
                 "dTonLdECMAHMSPowAvgBin",
                 "dTonLdECMAHMS05ExBin",
                 "dRoughECMAHMS10ExBin",
                 "dRoughECMAHMS05ExBin",
                 "dFluctHMS10ExBin",
                 "dFluctHMS05ExBin",
                 "dSharpAuresISO3PowAvgMaxLR",
                 "dSharpAuresISO305ExMaxLR",
                 "dImpulsHMSPowAvgMaxLR",
                 "dImpulsHMSAvgMaxLR",
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
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

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
    
    # fix problem filenames
    if stimulus.find("YnTy") != -1:
        stimulus = stimulus.replace("YnTy", "H520")
    
    if stimulus.find(".wav_Pa.wav") != -1:
        stimulus = stimulus.replace(".wav_Pa.wav", "_Pa.wav")
    
    # Calculate ECMA-418-2:2022 Sottek Hearing Model overall loudness from
    # 2-channel specific loudness
    # left channel
    specLoudHMSL = pd.DataFrame(workbookdata['Sheet1'].iloc[14:, 1:].values,
                                columns=workbookdata['Sheet1'].iloc[13, 1:],
                                index=workbookdata['Sheet1'].iloc[14:, 0])
    # right channel
    specLoudHMSR = pd.DataFrame(workbookdata['Sheet2'].iloc[14:, 1:].values,
                                columns=workbookdata['Sheet2'].iloc[13, 1:],
                                index=workbookdata['Sheet2'].iloc[14:, 0])
    # binaural specific loudness (ECMA-418-2:2022 Equation 118)
    specLoudHMSBin = ((specLoudHMSL**2
                        + specLoudHMSR**2)/2).pow(0.5)
    # binaural time-dependent loudness (ECMA-418-2:2022 Equation 116)
    loudHMSTimeVar = specLoudHMSBin.sum(axis=0)*0.5

    # mask for start/end skip
    loudHMSTimeVarMask = loudHMSTimeVar.loc[(loudHMSTimeVar.index.values
                                              > start_skipT)
                                            & (loudHMSTimeVar.index.values
                                                < loudHMSTimeVar.index.values.max()
                                                - end_skipT)]

    # binaural overall (power-averaged) loudness (ECMA-418-2:2022 Equation 117)
    loudHMSPowAvg = ((loudHMSTimeVarMask**(1/np.log10(2))).sum()
                      / len(loudHMSTimeVarMask))**np.log10(2)

    # Calculate ECMA-418-2:2022 Sottek Hearing Model overall tonality from
    # 2-channel specific tonality
    # left channel
    specTonalHMSL = pd.DataFrame(workbookdata['Sheet3'].iloc[14:, 1:].values,
                                  columns=workbookdata['Sheet3'].iloc[13, 1:],
                                  index=workbookdata['Sheet3'].iloc[14:, 0])
    # right channel
    specTonalHMSR = pd.DataFrame(workbookdata['Sheet4'].iloc[14:, 1:].values,
                                  columns=workbookdata['Sheet4'].iloc[13, 1:],
                                  index=workbookdata['Sheet4'].iloc[14:, 0])

    # 2-channel time-varing tonality (max, not integration)
    tonalHMSTimeVar = pd.concat([specTonalHMSL.max(axis=0),
                                 specTonalHMSR.max(axis=0)],
                                axis=1)
    
    # 2-channel time-varing tonality (integrated with adjustment to match 40 dB 1 kHz sine to 1 tu)
    tonalIntTimeVar = pd.concat([specTonalHMSL.sum(axis=0)*0.5*0.348088948583815,
                                 specTonalHMSR.sum(axis=0)*0.5*0.348088948583815],
                                axis=1)

    # mask for start/end skip and values <= 0.02
    tonalHMSTimeVarMaskL = tonalHMSTimeVar.loc[(tonalHMSTimeVar.index.values
                                                > start_skipT)
                                                & (tonalHMSTimeVar.index.values
                                                  < tonalHMSTimeVar.index.values.max()
                                                  - end_skipT)
                                                & (tonalHMSTimeVar.loc[:, 0].values > 0.02),
                                                0]
    tonalHMSTimeVarMaskR = tonalHMSTimeVar.loc[(tonalHMSTimeVar.index.values
                                                > start_skipT)
                                                & (tonalHMSTimeVar.index.values
                                                  < tonalHMSTimeVar.index.values.max()
                                                  - end_skipT)
                                                & (tonalHMSTimeVar.loc[:, 1].values > 0.02),
                                                1]
    
    # mask for start/end skip and values <= 0.02
    # NOTE: uses mask from HMS tonality
    tonalIntTimeVarMaskL = tonalIntTimeVar.loc[(tonalIntTimeVar.index.values
                                                > start_skipT)
                                                & (tonalIntTimeVar.index.values
                                                  < tonalIntTimeVar.index.values.max()
                                                  - end_skipT)
                                                & (tonalHMSTimeVar.loc[:, 0].values > 0.02),
                                                0]  # see NOTE above
    tonalIntTimeVarMaskR = tonalIntTimeVar.loc[(tonalIntTimeVar.index.values
                                                > start_skipT)
                                                & (tonalIntTimeVar.index.values
                                                  < tonalIntTimeVar.index.values.max()
                                                  - end_skipT)
                                                & (tonalHMSTimeVar.loc[:, 1].values > 0.02),
                                                1]  # see NOTE above
    
    # 2-channel time-averaged tonality (omitting T<=0.02)
    tonalHMSAvgL = tonalHMSTimeVarMaskL.mean(axis=0)
    tonalHMSAvgR = tonalHMSTimeVarMaskR.mean(axis=0)
    # max of L/R
    tonalHMSAvgMaxLR = max(tonalHMSAvgL, tonalHMSAvgR)
    
    # 2-channel 5% exceeded tonality (omitting T<=0.02)
    tonalHMS05ExL = tonalHMSTimeVarMaskL.quantile(q=0.95)
    tonalHMS05ExR = tonalHMSTimeVarMaskR.quantile(q=0.95)
    # max of L/R
    tonalHMS05ExMaxLR = max(tonalHMS05ExL, tonalHMS05ExR)
    
    # 2-channel time-averaged integated tonality (omitting T<=0.02)
    # NOTE: uses mask from HMS tonality
    tonalIntAvgL = tonalIntTimeVarMaskL.mean(axis=0)
    tonalIntAvgR = tonalIntTimeVarMaskR.mean(axis=0)
    # max of L/R
    tonalIntAvgMaxLR = max(tonalIntAvgL, tonalIntAvgR)
    
    # 2-channel 5% exceeded integrated tonality (omitting T<=0.02)
    tonalInt05ExL = tonalIntTimeVarMaskL.quantile(q=0.95)
    tonalInt05ExR = tonalIntTimeVarMaskR.quantile(q=0.95)
    # max of L/R
    tonalInt05ExMaxLR = max(tonalInt05ExL, tonalInt05ExR)

    # Calculate ECMA-418-2:2022 Sottek Hearing Model binaural overall roughness
    # from binaural specific roughness
    specRoughHMSBin = pd.DataFrame(workbookdata['Sheet5'].iloc[14:, 1:].values,
                                    columns=workbookdata['Sheet5'].iloc[13, 1:],
                                    index=workbookdata['Sheet5'].iloc[14:, 0])
    
    # binaural time-varying roughness
    roughHMSTimeVar = specRoughHMSBin.sum(axis=0)*0.5
    
    # mask for start/end skip
    roughHMSTimeVarMask = roughHMSTimeVar.loc[(roughHMSTimeVar.index.values
                                                > start_skipT)
                                              & (roughHMSTimeVar.index.values
                                                  < roughHMSTimeVar.index.values.max()
                                                  - end_skipT)]
    
    # binaural overall (90th percentile = 10% exceeded) roughness
    roughHMS10Ex = roughHMSTimeVarMask.quantile(q=0.9)
    # binaural overall (95th percentile = 5% exceeded) roughness
    roughHMS05Ex = roughHMSTimeVarMask.quantile(q=0.95)

    # Calculate overall Sottek Hearing Model fluctuation strength from
    # 2-channel specific fluctuation strength 
    # calculate according to ECMA-418-2:2022 approach for roughness
    # left channel
    specFluctHMSL = pd.DataFrame(workbookdata['Sheet6'].iloc[14:, 1:].values,
                                  columns=workbookdata['Sheet6'].iloc[13, 1:],
                                  index=workbookdata['Sheet6'].iloc[14:, 0])
    # right channel
    specFluctHMSR = pd.DataFrame(workbookdata['Sheet7'].iloc[14:, 1:].values,
                                  columns=workbookdata['Sheet7'].iloc[13, 1:],
                                  index=workbookdata['Sheet7'].iloc[14:, 0])
    
    # binaural specific fluctuation strength
    # (using ECMA-418-2:2022 Equation 112 for roughness)
    specFluctHMSBin = ((specFluctHMSL**2
                        + specFluctHMSR**2)/2).pow(0.5)
    # binaural time-dependent fluctuation strength
    # (using ECMA-418-2:2022 Equation 111 for roughness)
    fluctHMSTimeVar = specFluctHMSBin.sum(axis=0)*0.5
    
    # mask for start/end skip
    fluctHMSTimeVarMask = fluctHMSTimeVar.loc[(fluctHMSTimeVar.index.values
                                                > start_skipT)
                                              & (fluctHMSTimeVar.index.values
                                                  < fluctHMSTimeVar.index.values.max()
                                                  - end_skipT)]
    
    # binaural overall (90th percentile = 10% exceeded) fluctuation strength
    # (using ECMA-418-2:2022 Section 7.1.8 for roughness)
    fluctHMSTime10Ex = fluctHMSTimeVarMask.quantile(q=0.9)
    # binaural overall (95th percentile = 5% exceeded) fluctuation strength
    fluctHMSTime05Ex = fluctHMSTimeVarMask.quantile(q=0.95)

    # Calculate overall ISO 532-1 loudness from 2-channel time-varing loudness
    loudISO1TimeVar = pd.DataFrame(workbookdata['Sheet8'].iloc[13:, 1:3].values,
                                    columns=workbookdata['Sheet8'].iloc[12, 1:3],
                                    index=workbookdata['Sheet8'].iloc[13:, 0])
    
    # mask for start/end skip and 0 values
    loudISO1TimeVarMask = loudISO1TimeVar.loc[np.tile((loudISO1TimeVar.index.values
                                                      > start_skipT), (2, 1)).transpose()
                                              & np.tile((loudISO1TimeVar.index.values
                                                          < loudISO1TimeVar.index.values.max()
                                                          - end_skipT), (2, 1)).transpose()]
    
    # 2-channel overall (5% exceeded = 95th percentile) loudness
    loudISO105Ex = loudISO1TimeVarMask.quantile(q=0.95)
    # max of l/r channel (5% exceeded = 95th percentile) loudness
    loudISO105ExMaxLR = loudISO105Ex.max()

    # 2-channel overall (power-averaged) loudness
    loudISO1PowAvg = ((loudISO1TimeVarMask**(1/np.log10(2))).sum(axis=0)
                      / len(loudISO1TimeVarMask))**np.log10(2)
    # max of l/r channel (95th-percentile) loudness
    loudISO1PowAvgMaxLR = loudISO1PowAvg.max()

    # Calculate overall ISO 532-3 loudness from binaural time-varing loudness
    loudISO3TimeVar = pd.DataFrame(workbookdata['Sheet9'].iloc[13:, 1:2].values,
                                    columns=workbookdata['Sheet9'].iloc[12, 1:2],
                                    index=workbookdata['Sheet9'].iloc[13:, 0])

    # mask for start/end skip
    loudISO3TimeVarMask = loudISO3TimeVar.loc[(loudISO3TimeVar.index.values
                                                > start_skipT)
                                              & (loudISO3TimeVar.index.values
                                                  < loudISO3TimeVar.index.values.max()
                                                  - end_skipT)]

    # binaural overall (power-averaged) loudness
    loudISO3PowAvg = ((loudISO3TimeVarMask**(1/np.log10(2))).sum(axis=0)
                      / len(loudISO3TimeVarMask))**np.log10(2)
    loudISO3PowAvg = loudISO3PowAvg.iloc[0]  # convert series to float

    # Calculate overall Aures+ISO532-3 sharpness from 2-channel time-varying
    # sharpness
    sharpAISO3TimeVar = pd.DataFrame(workbookdata['Sheet10'].iloc[13:, 1:3].values,
                                      columns=workbookdata['Sheet10'].iloc[12, 1:3],
                                      index=workbookdata['Sheet10'].iloc[13:, 0])

    # mask for start/end skip
    sharpAISO3TimeVarMask = sharpAISO3TimeVar.loc[np.tile((sharpAISO3TimeVar.index.values
                                                            > start_skipT), (2, 1)).transpose()
                                                  & np.tile((sharpAISO3TimeVar.index.values
                                                              < sharpAISO3TimeVar.index.values.max()
                                                              - end_skipT), (2, 1)).transpose()]

    # 2-channel overall (power-averaged) sharpness
    sharpAISO3PowAvg = ((sharpAISO3TimeVarMask**(1/np.log10(2))).sum(axis=0)
                        / len(sharpAISO3TimeVarMask))**np.log10(2)
    # max of l/r channel overall (power-averaged) sharpness
    sharpAISO3PowAvgMaxLR = sharpAISO3PowAvg.max()
    # 2-channel overall 5% exceeded sharpness
    sharpAISO305Ex = sharpAISO3TimeVarMask.quantile(q=0.95)
    # max of l/r channel overall 5% exceeded sharpness
    sharpAISO305ExMaxLR = sharpAISO305Ex.max()
    
    # Calculate overall Sottek Hearing Model impulsiveness from 2-channel
    # time-varying impulsiveness
    impulsHMSTimeVar = pd.DataFrame(workbookdata['Sheet12'].iloc[13:, 1:3].values,
                                        columns=workbookdata['Sheet12'].iloc[12, 1:3],
                                        index=workbookdata['Sheet12'].iloc[13:, 0])

    # mask for start/end skip and 0 values
    impulsHMSTimeVarMask = impulsHMSTimeVar.loc[np.tile((impulsHMSTimeVar.index.values
                                                          > start_skipT), (2, 1)).transpose()
                                                & np.tile((impulsHMSTimeVar.index.values
                                                            < impulsHMSTimeVar.index.values.max()
                                                            - end_skipT), (2, 1)).transpose()]

    # 2-channel overall (power-averaged) impulsiveness
    impulsHMSPowAvg = ((impulsHMSTimeVarMask**(1/np.log10(2))).sum(axis=0)
                        / len(impulsHMSTimeVarMask))**np.log10(2)
    # max of l/r overall (power-averaged) impulsiveness
    impulsHMSPowAvgMaxLR = impulsHMSPowAvg.max()
    # 2-channel overall (mean) impulsiveness
    impulsHMSMean = impulsHMSTimeVarMask.mean(axis=0)
    # max of l/r channel overall (mean) impulsiveness
    impulsHMSMeanMaxLR = impulsHMSMean.max()

    # 2-channel overall 5% exceeded impulsiveness
    impulsHMS05Ex = impulsHMSTimeVarMask.quantile(q=0.95)

    # max of l/r channel overall 5% exceeded impulsiveness
    impulsHMS05ExMaxLR = impulsHMS05Ex.max()

    # Calculate ECMA-418-2:2022 Sottek Hearing Model overall tonal loudness from
    # 2-channel specific tonal loudness
    # left channel
    specTonLdHMSL = pd.DataFrame(workbookdata['Sheet13'].iloc[14:, 1:].values,
                                  columns=workbookdata['Sheet13'].iloc[13, 1:],
                                  index=workbookdata['Sheet13'].iloc[14:, 0])
    # right channel
    specTonLdHMSR = pd.DataFrame(workbookdata['Sheet14'].iloc[14:, 1:].values,
                                  columns=workbookdata['Sheet14'].iloc[13, 1:],
                                  index=workbookdata['Sheet14'].iloc[14:, 0])
    # binaural specific tonal loudness (ECMA-418-2:2022 Equation 118)
    specTonLdHMSBin = ((specTonLdHMSL**2
                        + specTonLdHMSR**2)/2).pow(0.5)
    # binaural time-dependent tonal loudness (ECMA-418-2:2022 Equation 116)
    tonLdHMSTimeVar = specTonLdHMSBin.sum(axis=0)*0.5

    # mask for start/end skip
    tonLdHMSTimeVarMask = tonLdHMSTimeVar.loc[(tonLdHMSTimeVar.index.values
                                              > start_skipT)
                                              & (tonLdHMSTimeVar.index.values
                                                  < tonLdHMSTimeVar.index.values.max()
                                                  - end_skipT)]

    # binaural overall (power-averaged) tonal loudness (ECMA-418-2:2022
    # Equation 117)
    tonLdHMSPowAvg = ((tonLdHMSTimeVarMask**(1/np.log10(2))).sum()
                      / len(tonLdHMSTimeVarMask))**np.log10(2)

    # binaural 5% exceeded tonal loudness
    tonLdHMS05Ex = tonLdHMSTimeVarMask.quantile(q=0.95)

    # add results to output DataFrame
    dataByStim.loc[stimulus, 'LoudECMAHMSPowAvgBin'] = loudHMSPowAvg
    dataByStim.loc[stimulus, 'TonalECMAHMSAvgMaxLR'] = tonalHMSAvgMaxLR
    dataByStim.loc[stimulus, 'TonalECMAHMS05ExMaxLR'] = tonalHMS05ExMaxLR
    dataByStim.loc[stimulus, 'TonalIntAvgMaxLR'] = tonalIntAvgMaxLR
    dataByStim.loc[stimulus, 'TonalInt05ExMaxLR'] = tonalInt05ExMaxLR
    dataByStim.loc[stimulus, 'TonLdECMAHMSPowAvgBin'] = tonLdHMSPowAvg
    dataByStim.loc[stimulus, 'TonLdECMAHMS05ExBin'] = tonLdHMS05Ex
    dataByStim.loc[stimulus, 'RoughECMAHMS10ExBin'] = roughHMS10Ex
    dataByStim.loc[stimulus, 'RoughECMAHMS05ExBin'] = roughHMS05Ex
    dataByStim.loc[stimulus, 'FluctHMS10ExBin'] = fluctHMSTime10Ex
    dataByStim.loc[stimulus, 'FluctHMS05ExBin'] = fluctHMSTime05Ex
    dataByStim.loc[stimulus, 'LoudISO105ExMaxLR'] = loudISO105ExMaxLR
    dataByStim.loc[stimulus, 'LoudISO1PowAvgMaxLR'] = loudISO1PowAvgMaxLR
    dataByStim.loc[stimulus, 'LoudISO3PowAvgBin'] = loudISO3PowAvg
    dataByStim.loc[stimulus, 'SharpAuresISO3PowAvgMaxLR'] = sharpAISO3PowAvgMaxLR
    dataByStim.loc[stimulus, 'SharpAuresISO305ExMaxLR'] = sharpAISO305ExMaxLR
    dataByStim.loc[stimulus, 'ImpulsHMSPowAvgMaxLR'] = impulsHMSPowAvgMaxLR
    dataByStim.loc[stimulus, 'ImpulsHMSAvgMaxLR'] = impulsHMSMeanMaxLR
    dataByStim.loc[stimulus, 'ImpulsHMS05ExMaxLR'] = impulsHMS05ExMaxLR

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
        ambSpecTonalHMSLMovAvg = specTonalHMSL.T.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean().T
        ambSpecTonalHMSRMovAvg = specTonalHMSR.T.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean().T
        ambSpecRoughHMSBinMovAvg = specRoughHMSBin.T.rolling(window=int(np.ceil(sampleRateRoughECMA*windowT))).mean().T
        ambSpecFluctHMSLMovAvg = specFluctHMSL.T.rolling(window=int(np.ceil((sampleRateFluctHMS)*windowT))).mean().T
        ambSpecFluctHMSRMovAvg = specFluctHMSR.T.rolling(window=int(np.ceil((sampleRateFluctHMS)*windowT))).mean().T
        ambSharpAISO3TimeVarMovAvg = sharpAISO3TimeVar.rolling(window=int(np.ceil(sampleRateLoudISO3*windowT))).mean()
        ambImpulsHMSTimeVarMovAvg = impulsHMSTimeVar.rolling(window=int(np.ceil(sampleRateImpulsHMS*windowT))).mean()
        ambSpecTonLdHMSLMovAvg = specTonLdHMSL.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()
        ambSpecTonLdHMSRMovAvg = specTonLdHMSR.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()

    elif stimulus[0:3] in ["A1_", "A2_", "B2_"]:
        # calculate moving average values for combined stimulus
        specTonalHMSLMovAvg = specTonalHMSL.T.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean().T
        specTonalHMSRMovAvg = specTonalHMSR.T.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean().T
        specRoughHMSBinMovAvg = specRoughHMSBin.T.rolling(window=int(np.ceil(sampleRateRoughECMA*windowT))).mean().T
        specFluctHMSLMovAvg = specFluctHMSL.T.rolling(window=int(np.ceil((sampleRateFluctHMS)*windowT))).mean().T
        specFluctHMSRMovAvg = specFluctHMSR.T.rolling(window=int(np.ceil((sampleRateFluctHMS)*windowT))).mean().T
        sharpAISO3TimeVarMovAvg = sharpAISO3TimeVar.rolling(window=int(np.ceil(sampleRateLoudISO3*windowT))).mean()
        impulsHMSTimeVarMovAvg = impulsHMSTimeVar.rolling(window=int(np.ceil(sampleRateImpulsHMS*windowT))).mean()
        specTonLdHMSLMovAvg = specTonLdHMSL.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()
        specTonLdHMRSMovAvg = specTonLdHMSR.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()

        # # calculate differences and make negative values 0
        dSpecTonalHMSL = np.maximum(specTonalHMSLMovAvg
                                    - ambSpecTonalHMSLMovAvg, 0)
        dSpecTonalHMSR = np.maximum(specTonalHMSRMovAvg
                                    - ambSpecTonalHMSRMovAvg, 0)
        dSpecRoughHMSBin = np.maximum(specRoughHMSBinMovAvg
                                      - ambSpecRoughHMSBinMovAvg, 0)
        dSpecFluctHMSL = np.maximum(specFluctHMSLMovAvg
                                    - ambSpecFluctHMSLMovAvg, 0)
        dSpecFluctHMSR = np.maximum(specFluctHMSRMovAvg
                                    - ambSpecFluctHMSRMovAvg, 0)
        dSharpAISO3TimeVar = np.maximum(sharpAISO3TimeVarMovAvg
                                        - ambSharpAISO3TimeVarMovAvg, 0)
        dImpulsHMSTimeVar = np.maximum(impulsHMSTimeVarMovAvg
                                        - ambImpulsHMSTimeVarMovAvg,
                                        0)
        dSpecTonLdHMSL = np.maximum(specTonLdHMSLMovAvg
                                    - ambSpecTonLdHMSLMovAvg, 0)
        dSpecTonLdHMSR = np.maximum(specTonLdHMRSMovAvg
                                    - ambSpecTonLdHMSRMovAvg, 0)

        # calculate aggregated difference values

        # 2-channel time-varing tonality
        dTonalHMSTimeVar = pd.concat([dSpecTonalHMSL.max(axis=0),
                                      dSpecTonalHMSR.max(axis=0)],
                                      axis=1)
        
        # 2-channel time-varing integrated tonality
        dTonalIntTimeVar = pd.concat([dSpecTonalHMSL.sum(axis=0)*0.5*0.348088948583815,
                                      dSpecTonalHMSR.sum(axis=0)*0.5*0.348088948583815],
                                      axis=1)

        # mask for start/end skip and values <= 0.02
        dTonalHMSTimeVarMaskL = dTonalHMSTimeVar.loc[(dTonalHMSTimeVar.index.values
                                                      > start_skipT)
                                                      & (dTonalHMSTimeVar.index.values
                                                        < dTonalHMSTimeVar.index.values.max()
                                                        - end_skipT)
                                                      & (dTonalHMSTimeVar.loc[:, 0].values > 0.02),
                                                      0]
        dTonalHMSTimeVarMaskR = dTonalHMSTimeVar.loc[(dTonalHMSTimeVar.index.values
                                                      > start_skipT)
                                                      & (dTonalHMSTimeVar.index.values
                                                        < dTonalHMSTimeVar.index.values.max()
                                                        - end_skipT)
                                                      & (dTonalHMSTimeVar.loc[:, 1].values > 0.02),
                                                      1]
        
        # mask for start/end skip and values <= 0.02
        # NOTE: uses mask from HMS tonality
        dTonalIntTimeVarMaskL = dTonalIntTimeVar.loc[(dTonalIntTimeVar.index.values
                                                      > start_skipT)
                                                      & (dTonalIntTimeVar.index.values
                                                        < dTonalIntTimeVar.index.values.max()
                                                        - end_skipT)
                                                      & (dTonalHMSTimeVar.loc[:, 0].values > 0.02),
                                                      0]  # see NOTE above
        dTonalIntTimeVarMaskR = dTonalIntTimeVar.loc[(dTonalIntTimeVar.index.values
                                                      > start_skipT)
                                                      & (dTonalIntTimeVar.index.values
                                                        < dTonalIntTimeVar.index.values.max()
                                                        - end_skipT)
                                                      & (dTonalHMSTimeVar.loc[:, 1].values > 0.02),
                                                      1]  # see NOTE above
        
        # 2-channel time-averaged tonality (omitting T<=0.02)
        dTonalHMSAvgL = dTonalHMSTimeVarMaskL.mean(axis=0)
        dTonalHMSAvgR = dTonalHMSTimeVarMaskL.mean(axis=0)
        dTonalHMSAvgMaxLR = max(dTonalHMSAvgL, dTonalHMSAvgR)
        # 2-channel 5% exceeded tonality (automatically omitting T<=0.02)
        dTonalHMS05ExL = dTonalHMSTimeVarMaskL.quantile(q=0.95)
        dTonalHMS05ExR = dTonalHMSTimeVarMaskL.quantile(q=0.95)
        dTonalHMS05ExMaxLR = max(dTonalHMS05ExL, dTonalHMS05ExR)

        # 2-channel time-averaged integrated tonality (omitting T<=0.02)
        dTonalIntAvgL = dTonalIntTimeVarMaskL.mean(axis=0)
        dTonalIntAvgR = dTonalIntTimeVarMaskL.mean(axis=0)
        dTonalIntAvgMaxLR = max(dTonalIntAvgL, dTonalIntAvgR)
        # 2-channel 5% exceeded integated tonality (automatically omitting T<=0.02)
        dTonalInt05ExL = dTonalIntTimeVarMaskL.quantile(q=0.95)
        dTonalInt05ExR = dTonalIntTimeVarMaskL.quantile(q=0.95)
        dTonalInt05ExMaxLR = max(dTonalInt05ExL, dTonalInt05ExR)

        # binaural time-varying roughness
        dRoughHMSTimeVar = dSpecRoughHMSBin.sum(axis=0)*0.5
        
        # mask for start/end skip
        dRoughHMSTimeVarMask = dRoughHMSTimeVar.loc[(dRoughHMSTimeVar.index.values
                                                      > start_skipT)
                                                    & (dRoughHMSTimeVar.index.values
                                                        < dRoughHMSTimeVar.index.values.max()
                                                        - end_skipT)]

        # binaural overall (90th percentile = 10% exceeded) roughness
        dRoughHMS10Ex = dRoughHMSTimeVarMask.quantile(q=0.9)
        # binaural overall (95th percentile = 5% exceeded) roughness
        dRoughHMS05Ex = dRoughHMSTimeVarMask.quantile(q=0.95)

        # binaural specific fluctuation strength
        # (using ECMA-418-2:2022 Equation 112 for roughness)
        dSpecFluctHMSBin = ((dSpecFluctHMSL**2
                              + dSpecFluctHMSR**2)/2).pow(0.2)
        # binaural time-dependent fluctuation strength
        # (using ECMA-418-2:2022 Equation 111 for roughness)
        dFluctHMSTimeVar = dSpecFluctHMSBin.sum(axis=0)*0.5

        # mask for start/end skip
        dFluctHMSTimeVarMask = dFluctHMSTimeVar.loc[(dFluctHMSTimeVar.index.values
                                                      > start_skipT)
                                                    & (dFluctHMSTimeVar.index.values
                                                        < dFluctHMSTimeVar.index.values.max()
                                                        - end_skipT)]
        
        # binaural overall (90th percentile = 10% exceeded) fluctuation strength
        # (using ECMA-418-2:2022 Section 7.1.8 for roughness)
        dFluctHMSTime10Ex = dFluctHMSTimeVarMask.quantile(q=0.9)
        # binaural overall (95th percentile = 5% exceeded) fluctuation strength
        dFluctHMSTime05Ex = dFluctHMSTimeVarMask.quantile(q=0.95)

        # 2-channel sharpness masked for start/end skip
        dSharpAISO3TimeVarMask = dSharpAISO3TimeVar.loc[np.tile((dSharpAISO3TimeVar.index.values
                                                                  > start_skipT), (2, 1)).transpose()
                                                        & np.tile((dSharpAISO3TimeVar.index.values
                                                                    < dSharpAISO3TimeVar.index.values.max()
                                                                    - end_skipT), (2, 1)).transpose()]

        # 2-channel overall (power-averaged) sharpness
        dSharpAISO3PowAvg = ((dSharpAISO3TimeVarMask**(1/np.log10(2))).sum(axis=0)
                              / len(dSharpAISO3TimeVarMask))**np.log10(2)
        # max of l/r channel overall (power-averaged) sharpness
        dSharpAISO3PowAvgMaxLR = dSharpAISO3PowAvg.max()
        # 2-channel overall 5% exceeded sharpness
        dSharpAISO305Ex = dSharpAISO3TimeVarMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded sharpness
        dSharpAISO305ExMaxLR = dSharpAISO305Ex.max()

        # 2-channel impulsiveness masked for start/end skip
        dImpulsHMSTimeVarMask = dImpulsHMSTimeVar.loc[np.tile((dImpulsHMSTimeVar.index.values
                                                                > start_skipT), (2, 1)).transpose()
                                                      & np.tile((dImpulsHMSTimeVar.index.values
                                                                  < dImpulsHMSTimeVar.index.values.max()
                                                                  - end_skipT), (2, 1)).transpose()]

        # 2-channel overall (power-averaged) impulsiveness
        dImpulsHMSPowAvg = ((dImpulsHMSTimeVarMask**(1/np.log10(2))).sum(axis=0)
                            / len(dImpulsHMSTimeVarMask))**np.log10(2)
        # max of l/r overall (power-averaged) impulsiveness
        dImpulsHMSPowAvgMaxLR = dImpulsHMSPowAvg.max()
        # 2-channel overall (mean) impulsiveness
        dImpulsHMSMean = dImpulsHMSTimeVarMask.mean(axis=0)
        # max of l/r channel overall (mean) impulsiveness
        dImpulsHMSMeanMaxLR = dImpulsHMSMean.max()

        # 2-channel overall 5% exceeded impulsiveness
        dImpulsHMS05Ex = dImpulsHMSTimeVarMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded impulsiveness
        dImpulsHMS05ExMaxLR = dImpulsHMS05Ex.max()

        
        # binaural specific tonal loudness (ECMA-418-2:2022 Equation 118)
        dSpecTonLdHMSBin = ((dSpecTonLdHMSL**2
                            + dSpecTonLdHMSR**2)/2).pow(0.5)
        # binaural time-dependent tonal loudness (ECMA-418-2:2022 Equation 116)
        dTonLdHMSTimeVar = dSpecTonLdHMSBin.sum(axis=0)*0.5

        # mask for start/end skip
        dTonLdHMSTimeVarMask = dTonLdHMSTimeVar.loc[(dTonLdHMSTimeVar.index.values
                                                    > start_skipT)
                                                    & (dTonLdHMSTimeVar.index.values
                                                        < dTonLdHMSTimeVar.index.values.max()
                                                        - end_skipT)]

        # binaural overall (power-averaged) tonal loudness (ECMA-418-2:2022
        # Equation 117)
        dTonLdHMSPowAvg = ((dTonLdHMSTimeVarMask**(1/np.log10(2))).sum()
                            / len(dTonLdHMSTimeVarMask))**np.log10(2)

        # binaural 5% exceeded tonal loudness
        dTonLdHMS05Ex = dTonLdHMSTimeVarMask.quantile(q=0.95)

        # add results to output DataFrame
        dataByStim.loc[stimulus, 'dTonalECMAHMSAvgMaxLR'] = dTonalHMSAvgMaxLR
        dataByStim.loc[stimulus, 'dTonalECMAHMS05ExMaxLR'] = dTonalHMS05ExMaxLR
        dataByStim.loc[stimulus, 'dTonalIntAvgMaxLR'] = dTonalIntAvgMaxLR
        dataByStim.loc[stimulus, 'dTonalInt05ExMaxLR'] = dTonalInt05ExMaxLR
        dataByStim.loc[stimulus, 'dTonLdECMAHMSPowAvgBin'] = dTonLdHMSPowAvg
        dataByStim.loc[stimulus, 'dTonLdECMAHMS05ExBin'] = dTonLdHMS05Ex
        dataByStim.loc[stimulus, 'dRoughECMAHMS10ExBin'] = dRoughHMS10Ex
        dataByStim.loc[stimulus, 'dRoughECMAHMS05ExBin'] = dRoughHMS05Ex
        dataByStim.loc[stimulus, 'dFluctHMS10ExBin'] = dFluctHMSTime10Ex
        dataByStim.loc[stimulus, 'dFluctHMS05ExBin'] = dFluctHMSTime05Ex
        dataByStim.loc[stimulus, 'dSharpAuresISO3PowAvgMaxLR'] = dSharpAISO3PowAvgMaxLR
        dataByStim.loc[stimulus, 'dSharpAuresISO305ExMaxLR'] = dSharpAISO305ExMaxLR
        dataByStim.loc[stimulus, 'dImpulsHMSPowAvgMaxLR'] = dImpulsHMSPowAvgMaxLR
        dataByStim.loc[stimulus, 'dImpulsHMSAvgMaxLR'] = dImpulsHMSMeanMaxLR
        dataByStim.loc[stimulus, 'dImpulsHMS05ExMaxLR'] = dImpulsHMS05ExMaxLR

# end of for loop over psychoacoustic metrics xlsx files


# --------------------------
# MATLAB metric calculations
# --------------------------

# open csv file selection dialog and assign filepaths to list
# PROJECT NOTE: the PNL calculation files are in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Analysis/PNL?csf=1&web=1&e=ISpB2W
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance()

fileExts = "*.csv"
file = list(QFileDialog.getOpenFileName(filter=fileExts))[0]

metricsMLab = pd.read_csv(file, header=0, index_col=0)

dataByStim = dataByStim.merge(metricsMLab, how='outer',
                              left_index=True, right_index=True)

# reorganise SQMs in dataframe
cols2find = ['LoudECMAHMSPowAvgBin', 'FluctHMS10ExBin',
             'SharpAuresISO3PowAvgMaxLR', 'RoughECMAHMS10ExBin']

reps = [metricsMLab.columns.str.contains("PNL").sum(),
        metricsMLab.columns.str.contains("Fluct").sum(),
        metricsMLab.columns.str.contains("Rough").sum(),
        metricsMLab.columns.str.contains("Tonal").sum()]

repcol = list()
for ii, col in enumerate(cols2find):
    repcol += [col]*reps[ii]

for iiMetric, metric in enumerate(metricsMLab.columns):
    popcol = dataByStim.pop(metric)
    loccol = dataByStim.columns.get_loc(repcol[iiMetric])
    dataByStim.insert(loc=loccol, column=metric, value=popcol)


# add UAS-only and ambient-only data to combined stimuli
# ------------------------------------------------------
indicesDiffPsycho = ["dTonalECMAHMSAvgMaxLR",
                     "dTonalECMAHMS05ExMaxLR",
                     "dTonalIntAvgMaxLR",
                     "dTonalInt05ExMaxLR",
                     "dTonLdECMAHMSPowAvgBin",
                     "dTonLdECMAHMS05ExBin",
                     "dRoughECMAHMS10ExBin",
                     "dRoughECMAHMS05ExBin",
                     "dFluctHMS10ExBin",
                     "dFluctHMS05ExBin",
                     "dSharpAuresISO3PowAvgMaxLR",
                     "dSharpAuresISO305ExMaxLR",
                     "dImpulsHMSPowAvgMaxLR",
                     "dImpulsHMSAvgMaxLR",
                     "dImpulsHMSMaxMaxLR",
                     "dImpulsHMS05ExMaxLR"]

indicesPsycho = indicesPsycho + list(metricsMLab.columns)[3:]

indicesAbsPsycho = [index for index in indicesPsycho
                    if index not in indicesDiffPsycho]

indices = indicesAcoustic + list(metricsMLab.columns)[:3] + indicesAbsPsycho

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
dataByStim['LAeqLAF90diff'] = dataByStim['UASLAeqMaxLR'] - dataByStim['AmbLAF90ExMaxLR']
dataByStim['LASmaxLAF90diff'] = dataByStim['UASLASmaxMaxLR'] - dataByStim['AmbLAF90ExMaxLR']
dataByStim['LASmaxLAF50diff'] = dataByStim['UASLASmaxMaxLR'] - dataByStim['AmbLAF50ExMaxLR']
dataByStim['LASmaxLAeqdiff'] = dataByStim['UASLASmaxMaxLR'] - dataByStim['AmbLAeqMaxLR']
dataByStim['PNLMLAeqdiff'] = dataByStim['UASPNLmaxMaxLR'] - dataByStim['AmbLAeqMaxLR']
dataByStim['PNLTMLAeqdiff'] = dataByStim['UASPNLTmaxMaxLR'] - dataByStim['AmbLAeqMaxLR']
dataByStim['EPNLLAeqdiff'] = dataByStim['UASEPNLMaxLR'] - dataByStim['AmbLAeqMaxLR']
dataByStim['PNLMLAF90diff'] = dataByStim['UASPNLmaxMaxLR'] - dataByStim['AmbLAF90ExMaxLR']
dataByStim['PNLTMLAF90diff'] = dataByStim['UASPNLTmaxMaxLR'] - dataByStim['AmbLAF90ExMaxLR']
dataByStim['EPNLLAF90diff'] = dataByStim['UASEPNLMaxLR'] - dataByStim['AmbLAF90ExMaxLR']
dataByStim['PNLMLAF50diff'] = dataByStim['UASPNLmaxMaxLR'] - dataByStim['AmbLAF50ExMaxLR']
dataByStim['PNLTMLAF50diff'] = dataByStim['UASPNLTmaxMaxLR'] - dataByStim['AmbLAF50ExMaxLR']
dataByStim['EPNLLAF50diff'] = dataByStim['UASEPNLMaxLR'] - dataByStim['AmbLAF50ExMaxLR']


# ---------------------
# Partial loudness data
# ---------------------

# open csv file selection dialog and assign filepaths to list
# PROJECT NOTE: the calculated files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Analysis/deeuu_loudness/output?csf=1&web=1&e=ZvblMt
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

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
                              columns=['PartLoudMGSTPowAvg',
                                       'PartLoudMGST05Ex'])
partLoudnessLT = pd.DataFrame(index=filenamesLT,
                              columns=['PartLoudMGLTPowAvg',
                                       'PartLoudMGLT05Ex'])
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
        partLoudnessST.loc[filenames[ii], 'PartLoudMGSTPowAvg'] = partialLoudnessPowAvg
        partLoudnessST.loc[filenames[ii], 'PartLoudMGST05Ex'] = partialLoudness05Ex
    elif file.find("Long") != -1:
        partLoudnessLT.loc[filenames[ii], 'PartLoudMGLTPowAvg'] = partialLoudnessPowAvg
        partLoudnessLT.loc[filenames[ii], 'PartLoudMGLT05Ex'] = partialLoudness05Ex

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
dataByStim.loc[NoUASStims, ['PartLoudMGSTPowAvg', 'PartLoudMGST05Ex',
                            'PartLoudMGLTPowAvg', 'PartLoudMGLT05Ex'] ] = 0

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
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

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

# add change in response data
partAResponses['dValence'] = np.nan
partAResponses['dArousal'] = np.nan
partAResponses['dAnnoyance'] = np.nan
for response in ['Valence', 'Arousal', 'Annoyance']:
    for iD in partAResponses['ID#'].unique():
        partAResponses.loc[(partAResponses['ID#']
                            == iD)
                           & (partAResponses['Stim File'].str.contains("A1")),
                           "d" + response] = (partAResponses.loc[(partAResponses['ID#']
                                                                  == iD)
                                                                 &
                                                                 (partAResponses['Stim File'].str.contains("A1"))].loc[:,
                                                                                                                       response]
                                              - partAResponses.loc[(partAResponses['ID#']
                                                                    == iD)
                                                                   & (partAResponses['Stim File']
                                                                      == 'A1.wav'),
                                                                   response].values)
        partAResponses.loc[(partAResponses['ID#']
                            == iD)
                           & (partAResponses['Stim File'].str.contains("A2")),
                           "d" + response] = (partAResponses.loc[(partAResponses['ID#']
                                                                  == iD)
                                                                 &
                                                                 (partAResponses['Stim File'].str.contains("A2"))].loc[:,
                                                                                                                       response]
                                              - partAResponses.loc[(partAResponses['ID#']
                                                                    == iD)
                                                                   & (partAResponses['Stim File']
                                                                      == 'A2.wav'),
                                                                   response].values)

    popcol = partAResponses.pop("d" + response)
    partAResponses.insert(loc=partAResponses.columns.get_loc('Recording'),
                          column="d" + response, value=popcol)

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

    partAdValence = partAResponses.loc[partAResponses['Recording'] == file,
                                       ['ID#', 'dValence']]
    columns = ["dValence_" + str(ID) for ID in partAdValence['ID#']]
    partAdValence = pd.DataFrame(data=np.array(partAdValence['dValence']),
                                 index=columns, columns=[file]).transpose()

    partAArousal = partAResponses.loc[partAResponses['Recording'] == file,
                                      ['ID#', 'Arousal']]
    columns = ["Arousal_" + str(ID) for ID in partAArousal['ID#']]
    partAArousal = pd.DataFrame(data=np.array(partAArousal['Arousal']),
                                index=columns, columns=[file]).transpose()

    partAdArousal = partAResponses.loc[partAResponses['Recording'] == file,
                                       ['ID#', 'dArousal']]
    columns = ["dArousal_" + str(ID) for ID in partAdArousal['ID#']]
    partAdArousal = pd.DataFrame(data=np.array(partAdArousal['dArousal']),
                                 index=columns, columns=[file]).transpose()

    partAAnnoy = partAResponses.loc[partAResponses['Recording'] == file,
                                    ['ID#', 'Annoyance']]
    columns = ["Annoyance_" + str(ID) for ID in partAAnnoy['ID#']]
    partAAnnoy = pd.DataFrame(data=np.array(partAAnnoy['Annoyance']),
                              index=columns, columns=[file]).transpose()

    partAdAnnoy = partAResponses.loc[partAResponses['Recording'] == file,
                                     ['ID#', 'dAnnoyance']]
    columns = ["dAnnoyance_" + str(ID) for ID in partAdAnnoy['ID#']]
    partAdAnnoy = pd.DataFrame(data=np.array(partAdAnnoy['dAnnoyance']),
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
    
    dvalenceSWtest = pd.DataFrame(data=[np.array(stats.shapiro(partAdValence))],
                                  columns=['dValenceSWtestW', 'dValenceSWtestp'],
                                  index=[file])
    darousalSWtest = pd.DataFrame(data=[np.array(stats.shapiro(partAdArousal))],
                                  columns=['dArousalSWtestW', 'dArousalSWtestp'],
                                  index=[file])
    dannoySWtest = pd.DataFrame(data=[np.array(stats.shapiro(partAdAnnoy))],
                                columns=['dAnnoySWtestW', 'dAnnoySWtestp'],
                                index=[file])

    valenceAgg = pd.DataFrame(data=[[np.percentile(partAValence.values,
                                                   q=50, axis=1,
                                                   method='median_unbiased')[0],
                                    np.mean(partAValence.values, axis=1)[0]]],
                              columns=['ValenceMedian', 'ValenceMean'],
                              index=[file])

    dvalenceAgg = pd.DataFrame(data=[[np.percentile(partAdValence.values,
                                                    q=50, axis=1,
                                                    method='median_unbiased')[0],
                                     np.mean(partAdValence.values, axis=1)[0]]],
                               columns=['dValenceMedian', 'dValenceMean'],
                               index=[file])

    arousalAgg = pd.DataFrame(data=[[np.percentile(partAArousal.values,
                                                   q=50, axis=1,
                                                   method='median_unbiased')[0],
                                    np.mean(partAArousal.values, axis=1)[0]]],
                              columns=['ArousalMedian', 'ArousalMean'],
                              index=[file])

    darousalAgg = pd.DataFrame(data=[[np.percentile(partAdArousal.values,
                                                    q=50, axis=1,
                                                    method='median_unbiased')[0],
                                     np.mean(partAdArousal.values, axis=1)[0]]],
                               columns=['dArousalMedian', 'dArousalMean'],
                               index=[file])

    annoyAgg = pd.DataFrame(data=[[np.percentile(partAAnnoy.values,
                                                 q=50, axis=1,
                                                 method='median_unbiased')[0],
                                  np.mean(partAAnnoy.values, axis=1)[0]]],
                            columns=['AnnoyMedian', 'AnnoyMean'],
                            index=[file])
    
    dannoyAgg = pd.DataFrame(data=[[np.percentile(partAdAnnoy.values,
                                                  q=50, axis=1,
                                                  method='median_unbiased')[0],
                                   np.mean(partAdAnnoy.values, axis=1)[0]]],
                             columns=['dAnnoyMedian', 'dAnnoyMean'],
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
                            partAdValence, partAdArousal, partAdAnnoy,
                            partAHighAnnoy,
                            partANotice, valenceAgg, arousalAgg,
                            annoyAgg, dvalenceAgg, darousalAgg,
                            dannoyAgg, highAnnoyAgg, noticeAgg])
        partAstats = partAstats.join([valenceSWtest, arousalSWtest,
                                      annoySWtest, dvalenceSWtest,
                                      darousalSWtest, dannoySWtest, ])
    else:
        partAstats.loc[file, valenceSWtest.columns] = valenceSWtest.loc[file]
        partAstats.loc[file, arousalSWtest.columns] = arousalSWtest.loc[file]
        partAstats.loc[file, annoySWtest.columns] = annoySWtest.loc[file]
        partAstats.loc[file, dvalenceSWtest.columns] = dvalenceSWtest.loc[file]
        partAstats.loc[file, darousalSWtest.columns] = darousalSWtest.loc[file]
        partAstats.loc[file, dannoySWtest.columns] = dannoySWtest.loc[file]

        partA.loc[file, partAValence.columns] = partAValence.loc[file]
        partA.loc[file, partAArousal.columns] = partAArousal.loc[file]
        partA.loc[file, partAAnnoy.columns] = partAAnnoy.loc[file]
        partA.loc[file, partAdValence.columns] = partAdValence.loc[file]
        partA.loc[file, partAdArousal.columns] = partAdArousal.loc[file]
        partA.loc[file, partAdAnnoy.columns] = partAdAnnoy.loc[file]
        partA.loc[file, partAHighAnnoy.columns] = partAHighAnnoy.loc[file]
        partA.loc[file, partANotice.columns] = partANotice.loc[file]
        partA.loc[file, valenceAgg.columns] = valenceAgg.loc[file]
        partA.loc[file, arousalAgg.columns] = arousalAgg.loc[file]
        partA.loc[file, annoyAgg.columns] = annoyAgg.loc[file]
        partA.loc[file, dvalenceAgg.columns] = dvalenceAgg.loc[file]
        partA.loc[file, darousalAgg.columns] = darousalAgg.loc[file]
        partA.loc[file, dannoyAgg.columns] = dannoyAgg.loc[file]
        partA.loc[file, highAnnoyAgg.columns] = highAnnoyAgg.loc[file]
        partA.loc[file, noticeAgg.columns] = noticeAgg.loc[file]


# Part B
# ------

# read in data and add column indicating the stimulus recording file
partBResponses = pd.read_excel(io=filepath, sheet_name="PtBResponse",
                               header=0)
partBResponses['Recording'] = [StimFile.split('.')[0] + "_CALBIN_Pa.wav"
                               for StimFile in partBResponses['Stim File']]

# add change in response data
partBResponses['dValence'] = np.nan
partBResponses['dArousal'] = np.nan
partBResponses['dAnnoyance'] = np.nan
for response in ['Valence', 'Arousal', 'Annoyance']:
    for iD in partBResponses['ID#'].unique():        
        partBResponses.loc[(partBResponses['ID#']
                            == iD)
                           & (partBResponses['Stim File'].str.contains("B2")),
                           "d" + response] = (partBResponses.loc[(partBResponses['ID#']
                                                                  == iD)
                                                                 &
                                                                 (partBResponses['Stim File'].str.contains("B2"))].loc[:,
                                                                                                                       response]
                                              - partBResponses.loc[(partBResponses['ID#']
                                                                    == iD)
                                                                   & (partBResponses['Stim File']
                                                                      == 'B2.wav'),
                                                                   response].values)

    popcol = partBResponses.pop("d" + response)
    partBResponses.insert(loc=partBResponses.columns.get_loc('Recording'),
                          column="d" + response, value=popcol)


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

    partBdValence = partBResponses.loc[partBResponses['Recording'] == file,
                                       ['ID#', 'dValence']]
    columns = ["dValence_" + str(ID) for ID in partBdValence['ID#']]
    partBdValence = pd.DataFrame(data=np.array(partBdValence['dValence']),
                                 index=columns, columns=[file]).transpose()

    partBArousal = partBResponses.loc[partBResponses['Recording'] == file,
                                      ['ID#', 'Arousal']]
    columns = ["Arousal_" + str(ID) for ID in partBArousal['ID#']]
    partBArousal = pd.DataFrame(data=np.array(partBArousal['Arousal']),
                                index=columns, columns=[file]).transpose()

    partBdArousal = partBResponses.loc[partBResponses['Recording'] == file,
                                       ['ID#', 'dArousal']]
    columns = ["dArousal_" + str(ID) for ID in partBdArousal['ID#']]
    partBdArousal = pd.DataFrame(data=np.array(partBdArousal['dArousal']),
                                 index=columns, columns=[file]).transpose()

    partBAnnoy = partBResponses.loc[partBResponses['Recording'] == file,
                                    ['ID#', 'Annoyance']]

    columns = ["Annoyance_" + str(ID) for ID in partBAnnoy['ID#']]
    partBAnnoy = pd.DataFrame(data=np.array(partBAnnoy['Annoyance']),
                              index=columns, columns=[file]).transpose()

    partBdAnnoy = partBResponses.loc[partBResponses['Recording'] == file,
                                     ['ID#', 'dAnnoyance']]

    columns = ["dAnnoyance_" + str(ID) for ID in partBdAnnoy['ID#']]
    partBdAnnoy = pd.DataFrame(data=np.array(partBdAnnoy['dAnnoyance']),
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
    
    dvalenceSWtest = pd.DataFrame(data=[np.array(stats.shapiro(partBdValence))],
                                  columns=['dValenceSWtestW', 'dValenceSWtestp'],
                                  index=[file])
    darousalSWtest = pd.DataFrame(data=[np.array(stats.shapiro(partBdArousal))],
                                  columns=['dArousalSWtestW', 'dArousalSWtestp'],
                                  index=[file])
    dannoySWtest = pd.DataFrame(data=[np.array(stats.shapiro(partBdAnnoy))],
                                columns=['dAnnoySWtestW', 'dAnnoySWtestp'],
                                index=[file])

    valenceAgg = pd.DataFrame(data=[[np.percentile(partBValence.values,
                                                   q=50, axis=1,
                                                   method='median_unbiased')[0],
                                    np.mean(partBValence.values, axis=1)[0]]],
                              columns=['ValenceMedian', 'ValenceMean'],
                              index=[file])
    
    dvalenceAgg = pd.DataFrame(data=[[np.percentile(partBdValence.values,
                                                    q=50, axis=1,
                                                    method='median_unbiased')[0],
                                     np.mean(partBdValence.values, axis=1)[0]]],
                               columns=['dValenceMedian', 'dValenceMean'],
                               index=[file])

    arousalAgg = pd.DataFrame(data=[[np.percentile(partBArousal.values,
                                                   q=50, axis=1,
                                                   method='median_unbiased')[0],
                                    np.mean(partBArousal.values, axis=1)[0]]],
                              columns=['ArousalMedian', 'ArousalMean'],
                              index=[file])

    darousalAgg = pd.DataFrame(data=[[np.percentile(partBdArousal.values,
                                                    q=50, axis=1,
                                                    method='median_unbiased')[0],
                                     np.mean(partBdArousal.values, axis=1)[0]]],
                               columns=['dArousalMedian', 'dArousalMean'],
                               index=[file])

    annoyAgg = pd.DataFrame(data=[[np.percentile(partBAnnoy.values,
                                                 q=50, axis=1,
                                                 method='median_unbiased')[0],
                                  np.mean(partBAnnoy.values, axis=1)[0]]],
                            columns=['AnnoyMedian', 'AnnoyMean'],
                            index=[file])

    dannoyAgg = pd.DataFrame(data=[[np.percentile(partBdAnnoy.values,
                                                  q=50, axis=1,
                                                  method='median_unbiased')[0],
                                   np.mean(partBdAnnoy.values, axis=1)[0]]],
                             columns=['dAnnoyMedian', 'dAnnoyMean'],
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
                            partBdValence, partBdArousal, partBdAnnoy,
                            partBHighAnnoy,
                            valenceAgg, arousalAgg, annoyAgg,
                            dvalenceAgg, darousalAgg, dannoyAgg, highAnnoyAgg])
        partBstats = partBstats.join([valenceSWtest, arousalSWtest,
                                      annoySWtest, dvalenceSWtest,
                                      darousalSWtest, dannoySWtest])
    else:
        partBstats.loc[file, valenceSWtest.columns] = valenceSWtest.loc[file]
        partBstats.loc[file, arousalSWtest.columns] = arousalSWtest.loc[file]
        partBstats.loc[file, annoySWtest.columns] = annoySWtest.loc[file]
        partBstats.loc[file, dvalenceSWtest.columns] = dvalenceSWtest.loc[file]
        partBstats.loc[file, darousalSWtest.columns] = darousalSWtest.loc[file]
        partBstats.loc[file, dannoySWtest.columns] = dannoySWtest.loc[file]
        
        partB.loc[file, partBValence.columns] = partBValence.loc[file]
        partB.loc[file, partBArousal.columns] = partBArousal.loc[file]
        partB.loc[file, partBAnnoy.columns] = partBAnnoy.loc[file]
        partB.loc[file, partBdValence.columns] = partBdValence.loc[file]
        partB.loc[file, partBdArousal.columns] = partBdArousal.loc[file]
        partB.loc[file, partBdAnnoy.columns] = partBdAnnoy.loc[file]
        partB.loc[file, partBHighAnnoy.columns] = partBHighAnnoy.loc[file]
        partB.loc[file, valenceAgg.columns] = valenceAgg.loc[file]
        partB.loc[file, arousalAgg.columns] = arousalAgg.loc[file]
        partB.loc[file, annoyAgg.columns] = annoyAgg.loc[file]
        partB.loc[file, dvalenceAgg.columns] = dvalenceAgg.loc[file]
        partB.loc[file, darousalAgg.columns] = darousalAgg.loc[file]
        partB.loc[file, dannoyAgg.columns] = dannoyAgg.loc[file]
        partB.loc[file, highAnnoyAgg.columns] = highAnnoyAgg.loc[file]
    
allResponses = pd.concat([partB, partA], axis=0, join='outer')
allResponses.sort_index(inplace=True)

# merge Shapiro-Wilks test results
allShapWilksTest = pd.concat([partAstats, partBstats]).sort_index()

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

dataByStim.insert(loc=8, column='UASType',
                  value=dataByStim.pop('UASType'))

# ambient sound environment
dataByStim.loc[np.logical_or(dataByStim.index.str.find("A2") != -1,
                             dataByStim.index.str.find("B2") != -1),
               'AmbientEnv'] = "Park"
dataByStim.loc[dataByStim.index.str.find("A1") != -1, 'AmbientEnv'] = "Street"

dataByStim.insert(loc=9, column='AmbientEnv',
                  value=dataByStim.pop('AmbientEnv'))


# ----------------------------------
# Prepare outputs for saving to file
# ----------------------------------


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

# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

outFilePath = QFileDialog.getExistingDirectory(caption="Choose output folder to save processed files")

dataByStimTest.to_csv(os.path.join(outFilePath,
                                   "refmap_listest1_testdata_ByStim.csv"))

dataByStimTestA.to_csv(os.path.join(outFilePath,
                                    "refmap_listest1_testdataA_ByStim.csv"))

dataByStimTestB.to_csv(os.path.join(outFilePath,
                                    "refmap_listest1_testdataB_ByStim.csv"))

dataByStimAux.to_csv(os.path.join(outFilePath,
                                  "refmap_listest1_auxdata.csv"))

dataByStim.to_csv(os.path.join(outFilePath,
                               "refmap_listest1_alldata_ByStim.csv"))

# save statistical tests
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

statsoutFilePath = QFileDialog.getExistingDirectory(caption="Choose output folder to save statistical analysis files")
allShapWilksTest.to_csv(os.path.join(statsoutFilePath,
                                     "refmap_listest1_alltestShapWilks.csv"))

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

allDataBySubj.to_csv(os.path.join(outFilePath,
                                  "refmap_listest1_testdata_BySubj.csv"),
                     index=False)

partADataBySubj.to_csv(os.path.join(outFilePath,
                                    "refmap_listest1_testdataA_BySubj.csv"),
                       index=False)

partBDataBySubj.to_csv(os.path.join(outFilePath,
                                    "refmap_listest1_testdataB_BySubj.csv"),
                       index=False)

# Re-save filtered datasets for noticeability analysis only
omitParticipants = [2, 5, 21, 32, 34, 36, 38, 39, 40, 44, 45]
omitColumns = ["UAS_noticed_" + str(partID) for partID in omitParticipants]

partADataBySubjNotice = partADataBySubj.loc[~partADataBySubj['ID#'].isin(omitParticipants), :]

partADataBySubjNotice.to_csv(os.path.join(outFilePath,
                                          "refmap_listest1_testdataANoticeFilt_BySubj.csv"),
                             index=False)

omitColumns = omitColumns + (["Arousal_" + str(partID) for partID in omitParticipants])
omitColumns = omitColumns + (["Valence_" + str(partID) for partID in omitParticipants])
omitColumns = omitColumns + (["Annoyance_" + str(partID) for partID in omitParticipants])
omitColumns = omitColumns + (["dArousal_" + str(partID) for partID in omitParticipants])
omitColumns = omitColumns + (["dValence_" + str(partID) for partID in omitParticipants])
omitColumns = omitColumns + (["dAnnoyance_" + str(partID) for partID in omitParticipants])
omitColumns = omitColumns + (["HighAnnoy_" + str(partID) for partID in omitParticipants])
dataByStimTestANotice = dataByStimTestA.drop(labels=omitColumns, axis=1)
dataByStimTestANotice.drop(labels=['ArousalMean', 'ArousalMedian',
                                   'ValenceMean', 'ValenceMedian',
                                   'AnnoyMean', 'AnnoyMedian',
                                   'dArousalMean', 'dArousalMedian',
                                   'dValenceMean', 'dValenceMedian',
                                   'dAnnoyMean', 'dAnnoyMedian',
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
dataByStimTestANotice['ArousalMedianFilt'] = np.percentile(dataByStimTestANotice[keepColumns],
                                                           q=50, axis=1,
                                                           method='median_unbiased')
keepColumns = [label.replace("Arousal_", "Valence_") for label in keepColumns]
dataByStimTestANotice['ValenceMeanFilt'] = dataByStimTestANotice[keepColumns].mean(axis=1)
dataByStimTestANotice['ValenceMedianFilt'] = np.percentile(dataByStimTestANotice[keepColumns],
                                                           q=50, axis=1,
                                                           method='median_unbiased')
keepColumns = [label.replace("Valence_", "Annoyance_") for label in keepColumns]
dataByStimTestANotice['AnnoyMeanFilt'] = dataByStimTestANotice[keepColumns].mean(axis=1)
dataByStimTestANotice['AnnoyMedianFilt'] = np.percentile(dataByStimTestANotice[keepColumns],
                                                         q=50, axis=1,
                                                         method='median_unbiased')
keepColumns = [label.replace("Annoyance_", "dArousal_") for label in keepColumns]
dataByStimTestANotice['dArousalMeanFilt'] = dataByStimTestANotice[keepColumns].mean(axis=1)
dataByStimTestANotice['dArousalMedianFilt'] = np.percentile(dataByStimTestANotice[keepColumns],
                                                            q=50, axis=1,
                                                            method='median_unbiased')
keepColumns = [label.replace("dArousal_", "dValence_") for label in keepColumns]
dataByStimTestANotice['dValenceMeanFilt'] = dataByStimTestANotice[keepColumns].mean(axis=1)
dataByStimTestANotice['dValenceMedianFilt'] = np.percentile(dataByStimTestANotice[keepColumns],
                                                            q=50, axis=1,
                                                            method='median_unbiased')
keepColumns = [label.replace("dValence_", "dAnnoyance_") for label in keepColumns]
dataByStimTestANotice['dAnnoyMeanFilt'] = dataByStimTestANotice[keepColumns].mean(axis=1)
dataByStimTestANotice['dAnnoyMedianFilt'] = np.percentile(dataByStimTestANotice[keepColumns],
                                                          q=50, axis=1,
                                                          method='median_unbiased')
keepColumns = [label.replace("dAnnoyance_", "HighAnnoy_") for label in keepColumns]
dataByStimTestANotice['HighAnnoyTotalFilt'] = dataByStimTestANotice[keepColumns].sum(axis=1)
dataByStimTestANotice['HighAnnoyPropFilt'] = dataByStimTestANotice['HighAnnoyTotalFilt']/len(keepParticipants)

dataByStimTestANotice.to_csv(os.path.join(outFilePath,
                                          "refmap_listest1_testdataANoticeFilt_ByStim.csv"))

# save pre-test and post-test response data to separate files

preTestResponses.to_csv(os.path.join(outFilePath,
                                     "refmap_listest1_pretestdata.csv"),
                        index=False)

postTestResponses.to_csv(os.path.join(outFilePath,
                                      "refmap_listest1_posttestdata.csv"),
                         index=False)

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
partAdAnnoyDataBySubjWide = partADataBySubj.pivot(index='ID#',
                                                  columns='Recording',
                                                  values='dAnnoyance')
partAdValenceDataBySubjWide = partADataBySubj.pivot(index='ID#',
                                                    columns='Recording',
                                                    values='dValence')
partAdArousalDataBySubjWide = partADataBySubj.pivot(index='ID#',
                                                    columns='Recording',
                                                    values='dArousal')
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
partAdAnnoyDataBySubjWide = partAdAnnoyDataBySubjWide.merge(prePostTestResponses,
                                                          on='ID#')
partAdValenceDataBySubjWide = partAdValenceDataBySubjWide.merge(prePostTestResponses,
                                                              on='ID#')
partAdArousalDataBySubjWide = partAdArousalDataBySubjWide.merge(prePostTestResponses,
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
partAdAnnoyDataBySubjWide.to_csv(os.path.join(outFilePath,
                                              "refmap_listest1_dAnnoyDataABySubjWide.csv"),
                                              index=False)
partAdValenceDataBySubjWide.to_csv(os.path.join(outFilePath,
                                                "refmap_listest1_dValenceDataABySubjWide.csv"),
                                   index=False)
partAdArousalDataBySubjWide.to_csv(os.path.join(outFilePath,
                                                "refmap_listest1_dArousalDataABySubjWide.csv"),
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
partBdAnnoyDataBySubjWide = partBDataBySubj.pivot(index='ID#',
                                                  columns='Recording',
                                                  values='dAnnoyance')
partBdValenceDataBySubjWide = partBDataBySubj.pivot(index='ID#',
                                                    columns='Recording',
                                                    values='dValence')
partBdArousalDataBySubjWide = partBDataBySubj.pivot(index='ID#',
                                                    columns='Recording',
                                                    values='dArousal')

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
partBdAnnoyDataBySubjWide = partBdAnnoyDataBySubjWide.merge(prePostTestResponses,
                                                           on='ID#')
partBdValenceDataBySubjWide = partBdValenceDataBySubjWide.merge(prePostTestResponses,
                                                               on='ID#')
partBdArousalDataBySubjWide = partBdArousalDataBySubjWide.merge(prePostTestResponses,
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