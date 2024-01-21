# -*- coding: utf-8 -*-

# script

# import statements
import sys
import os
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import QFileDialog, QApplication


# open xlsx file selection dialog and assign filepaths to list
app = QApplication(sys.argv)
fileExts = "*.xlsx"
filelist = list(QFileDialog.getOpenFileNames(filter=fileExts))[0]

# setup output DataFrame
indices = ["ECMA_HMSloudnessBin",
           "ECMA_HMStonalityMaxLR",
           "ECMA_HMStonalityMeanLR",
           "ECMA_HMSroughnessBin",
           ]
# output


# skip signal start and end for time-aggregation
start_skipT = 0.5
end_skipT = 0.5
# loop over files to analyse
for file in filelist:

    stimulus = pd.read_excel(io=file, sheet_name="Sheet1", usecols='F',
                             nrows=1,
                             header=None).iloc[0, 0].split(sep='\'')[1]
    # Calculate ECMA-418-2:2022 Sottek Hearing Model overall loudness from
    # 2-channel specific loudness
    # left channel
    specificLoudnessHMSL = pd.read_excel(io=file, sheet_name="Sheet1",
                                         skiprows=14, header=0,
                                         index_col=0).transpose()
    # right channel
    specificLoudnessHMSR = pd.read_excel(io=file, sheet_name="Sheet2",
                                         skiprows=14, header=0,
                                         index_col=0).transpose()
    # binaural specific loudness (ECMA-418-2:2022 Equation 118)
    specificLoudnessHMSBin = ((specificLoudnessHMSL**2
                               + specificLoudnessHMSR**2)/2).pow(1/2)
    # binaural time-dependent loudness (ECMA-418-2:2022 Equation 116)
    loudnessHMSTimeVar = specificLoudnessHMSBin.sum(axis=1)*0.5
    # binaural overall (power-averaged) loudness (ECMA-418-2:2022 Equation 117)
    loudnessHMSPowAvg = ((loudnessHMSTimeVar.iloc[int(np.ceil(187.5*start_skipT)):
                                                  -int(np.ceil(187.5*end_skipT))]**(1/np.log10(2))).sum()
                         / len(loudnessHMSTimeVar.iloc[int(np.ceil(187.5*start_skipT)):
                                                       -int(np.ceil(187.5*end_skipT))]))**np.log10(2)

    # Calculate ECMA-418-2:2022 Sottek Hearing Model overall tonality from
    # 2-channel specific tonality
    # left channel
    specificTonalityHMSL = pd.read_excel(io=file, sheet_name="Sheet3",
                                         skiprows=14, header=0,
                                         index_col=0).transpose()
    # right channel
    specificTonalityHMSR = pd.read_excel(io=file, sheet_name="Sheet4",
                                         skiprows=14, header=0,
                                         index_col=0).transpose()
    # 2-channel time-varing tonality (max, not integration)
    tonalityHMSTimeVar = pd.concat([specificTonalityHMSL.max(axis=1),
                                    specificTonalityHMSR.max(axis=1)],
                                   axis=1)
    # 2-channel time-averaged tonality (omitting T<=0.02)
    tonalityHMSAvgL = tonalityHMSTimeVar.iloc[int(np.ceil(187.5*start_skipT)):,
                                              0][tonalityHMSTimeVar.iloc[int(np.ceil(187.5*start_skipT)):
                                                                         -int(np.ceil(187.5*end_skipT)), 0]
                                                 >= 0.02].mean()
    tonalityHMSAvgR = tonalityHMSTimeVar.iloc[int(np.ceil(187.5*start_skipT)):,
                                              1][tonalityHMSTimeVar.iloc[int(np.ceil(187.5*start_skipT)):
                                                                         -int(np.ceil(187.5*end_skipT)), 1]
                                                 >= 0.02].mean()
    tonalityHMSAvgMaxLR = max(tonalityHMSAvgL, tonalityHMSAvgR)

    # Calculate ECMA-418-2:2022 Sottek Hearing Model binaural overall roughness
    # from binaural specific roughness
    specificRoughnessHMSBin = pd.read_excel(io=file, sheet_name="Sheet5",
                                            skiprows=14, header=0,
                                            index_col=0).transpose()
    # binaural time-varying roughness
    roughnessHMSTimeVar = specificRoughnessHMSBin.sum(axis=1)*0.5
    # binaural 90th percentile roughness
    roughnessHMSPc90 = roughnessHMSTimeVar.iloc[int(np.ceil(50*start_skipT)):
                                                -int(np.ceil(50*end_skipT))].quantile(0.9)
  
    # Calculate overall Sottek Hearing Model fluctuation strength from
    # 2-channel specific fluctuation strength 
    # calculate according to ECMA-418-2:2022 approach for roughness
    # left channel
    specificFluctStrHMSL = pd.read_excel(io=file, sheet_name="Sheet6",
                                         skiprows=14, header=0,
                                         index_col=0).transpose()
    # right channel
    specificFluctStrHMSR = pd.read_excel(io=file, sheet_name="Sheet7",
                                         skiprows=14, header=0,
                                         index_col=0).transpose()
    # binaural specific fluctuation strength
    # (using ECMA-418-2:2022 Equation 112 for roughness)
    specificFluctStrHMSBin = ((specificFluctStrHMSL**2
                               + specificFluctStrHMSR**2)/2).pow(1/2)
    # binaural time-dependent fluctuation strength
    # (using ECMA-418-2:2022 Equation 111 for roughness)
    fluctStrHMSTimeVar = specificFluctStrHMSBin.sum(axis=1)*0.5
    # binaural overall (90th-percentile) fluctuation strength
    # (using ECMA-418-2:2022 Section 7.1.8 for roughness)
    fluctStrHMSTimePc90 = fluctStrHMSTimeVar.iloc[int(np.ceil(start_skipT*(229390681/4e6))):
                                                  -int(np.ceil(end_skipT*(229390681/4e6)))].quantile(0.9)

    # Calculate overall ISO 532-1 loudness from 2-channel time-varing loudness
    loudnessISO1TimeVar = pd.read_excel(io=file, sheet_name="Sheet8",
                                        usecols='A:C', skiprows=13, header=0,
                                        index_col=0)
    # 2-channel overall (95th-percentile) loudness
    loudnessISO1Pc95 = loudnessISO1TimeVar.iloc[int(start_skipT*1e3):-int(end_skipT*1e3)].quantile(0.95)
    # max of l/r channel (95th-percentile) loudness
    loudnessISO1Pc95MaxLR = loudnessISO1Pc95.max()

    # 2-channel overall (power-averaged) loudness
    loudnessISO1PowAvg = ((loudnessISO1TimeVar.iloc[int(start_skipT*500):
                                                    -int(end_skipT*500)]**(1/np.log10(2))).sum()
                          / len(loudnessISO1TimeVar.iloc[int(start_skipT*500):
                                                         -int(end_skipT*500)]))**np.log10(2)
    # max of l/r channel (95th-percentile) loudness
    loudnessISO1PowAvgMaxLR = loudnessISO1PowAvg.max()

    # Calculate overall ISO 532-3 loudness from binaural time-varing loudness
    loudnessISO3TimeVar = pd.read_excel(io=file, sheet_name="Sheet9",
                                        usecols='A:B', skiprows=13, header=0,
                                        index_col=0)
    # binaural overall (power-averaged) loudness
    loudnessISO3PowAvg = ((loudnessISO3TimeVar.iloc[int(start_skipT*1e3):
                                                    -int(end_skipT*1e3)]**(1/np.log10(2))).sum()
                          / len(loudnessISO3TimeVar.iloc[int(start_skipT*1e3):
                                                         -int(end_skipT*1e3)]))**np.log10(2)

    # Calculate overall Aures+ISO532-3 sharpness from 2-channel time-varying
    # sharpness
    sharpnessAISO3TimeVar = pd.read_excel(io=file, sheet_name="Sheet10",
                                          usecols='A:C', skiprows=13, header=0,
                                          index_col=0)
    # 2-channel overall (mean) sharpness
    sharpnessAISO3Mean = sharpnessAISO3TimeVar.mean()
    # max of l/r channel overall (mean) sharpness
    sharpnessAISO3MeanMaxLR = sharpnessAISO3Mean.max()
    # 2-channel overall (power-averaged) sharpness
    sharpnessAISO3PowAvg = ((sharpnessAISO3TimeVar.iloc[int(start_skipT*1e3):
                                                        -int(end_skipT*1e3)]**(1/np.log10(2))).sum()
                            / len(sharpnessAISO3TimeVar.iloc[int(start_skipT*1e3):
                                                            -int(end_skipT*1e3)]))**np.log10(2)
    # max of l/r channel overall (power-averaged) sharpness
    sharpnessAISO3PowAvgMaxLR = sharpnessAISO3PowAvg.max()
    
    # Calculate overall Sottek Hearing Model impulsiveness from 2-channel
    # time-varying impulsiveness
    impulsivenessHMSTimeVar = pd.read_excel(io=file, sheet_name="Sheet12",
                                              usecols='A:C', skiprows=13, header=0,
                                              index_col=0)
    # 2-channel overall (power-averaged) impulsiveness
    impulsivenessHMSPowAvg = ((impulsivenessHMSTimeVar.iloc[int(np.ceil(60000/55*start_skipT)):
                                                            -int(np.ceil(60000/55*end_skipT))]**(1/np.log10(2))).sum()
                              / len(impulsivenessHMSTimeVar.iloc[int(np.ceil(60000/55*start_skipT)):-int(np.ceil(60000/55*end_skipT))]))**np.log10(2)
    # max of l/r overall (power-averaged) impulsiveness
    impulsivenessHMSPowAvgMaxLR = impulsivenessHMSPowAvg.max()
    # 2-channel overall (mean) impulsiveness
    impulsivenessHMSMean = impulsivenessHMSTimeVar.mean()
    # max of l/r channel overall (mean) impulsiveness
    impulsivenessHMSMeanMaxLR = impulsivenessHMSMean.max()
    # max of l/r channel overall (max) impulsiveness
    impulsivenessHMSMaxMaxLR = impulsivenessHMSTimeVar.max(axis=None)


# end of for loop over xlsx files


# open ascii file selection dialog and assign filepaths to list
app = QApplication(sys.argv)
fileExts = "*.asc"
filelist = list(QFileDialog.getOpenFileNames(filter=fileExts))[0]
filelist.sort()

# first, extract specific loudness for both ambient scenes
specificLoudnessISO3