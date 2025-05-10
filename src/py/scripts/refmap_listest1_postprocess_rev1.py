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
from scipy import stats, io
from warnings import simplefilter

# suppress pandas performance warnings
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

# enable copy-on-write mode for Pandas (will be default from Pandas 3.0)
pd.options.mode.copy_on_write = True

# ---------------------------------------------------------
# Initialise data frame for metric calculations and results
# ---------------------------------------------------------

# skip signal start and end for time-aggregation
start_skipT = 0.5
end_skipT = 0.5

# open wav file selection dialog and assign filepaths to list
# PROJECT NOTE: the calibrated files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Calibration/Post_calib_recs/Pa_calib?csf=1&web=1&e=7eX89P
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance()

fileExts = "*.wav"
filelist = list(QFileDialog.getOpenFileNames(caption="Open recording files in '03 Experiment\Experiment 1\Calibration\Post_calib_recs\Pa_calib'",
                                             filter=fileExts))[0]
filelist.sort()
filenames = [filepath.split('/')[-1] for filepath in filelist]
stemNames = [filename.replace("_CALBIN_Pa.wav", "") for filename in filenames]
sqmNames = [filename.replace("CALBIN", "CALHEQ") for filename in filenames
            if filename.find("G57") == -1]
sqmNames = sqmNames + [""]

dataByStim = pd.DataFrame(index=stemNames, dtype=float)

dataByStim['CALBINRecFiles'] = filenames
dataByStim['CALHEQRecFiles'] = sqmNames

# -------------------------------------
# Add categorical variables for stimuli
# -------------------------------------

# NOTE: the order of this section matters, as some categorical variable
# definitions depend on Boolean logic derived from other categorical variables
# added

# Experiment session part A or B
dataByStim.loc[(dataByStim.index.str.find("A1") == 0)
               | (dataByStim.index.str.find("A2") == 0)
               | (dataByStim.index.str.find("A_") == 0),
               'SessionPart'] = "A"
dataByStim.loc[(dataByStim.index.str.find('B2') == 0)
               | (dataByStim.index.str.find("B_") == 0),
               'SessionPart'] = "B"

popcol = dataByStim.pop('SessionPart')
dataByStim.insert(loc=0, column='SessionPart',
                  value=popcol)

# Stimulus duration
for file in dataByStim.index:
    if file.split('_')[0].find("A") == 0:
        dataByStim.loc[file, 'StimDuration'] = 25
    elif file.split('_')[0].find("B") == 0:
        dataByStim.loc[file, 'StimDuration'] = 75
    else:
        dataByStim.loc[file, 'StimDuration'] = np.nan

# UAS LAeq level category
# Add the No UAS entries first
dataByStim.loc[(dataByStim.CALBINRecFiles.str.find("B2_CALBIN") != -1)
               | (dataByStim.CALBINRecFiles.str.find("A1_CALBIN") != -1)
               | (dataByStim.CALBINRecFiles.str.find("A2_CALBIN") != -1),
               'UASLAeq'] = "No UAS"

dataByStim.loc[(dataByStim.index.str.find("_F_1") != -1)
               | ((dataByStim.index.str.find("_1") != -1)
                  & (dataByStim.index.str.find("_F_2") == -1)),
               'UASLAeq'] = int(60)

dataByStim.loc[(dataByStim.index.str.find("_F_2") != -1)
               | ((dataByStim.index.str.find("_2") != -1)
                  & (dataByStim.index.str.find("_F_1") == -1)),
               'UASLAeq'] = int(54)

dataByStim.loc[(dataByStim.index.str.find("_3") != -1)
               & ~((dataByStim.index.str.find("_F_1") != -1)
                   | (dataByStim.index.str.find("_F_2") != -1)),
               'UASLAeq'] = int(48)
dataByStim.loc[dataByStim.index.str.find("_4") != -1, 'UASLAeq'] = int(42)

# ambient sound LAeq level category
dataByStim.loc[dataByStim.index.str.find("A1") != -1, 'AmbientLAeq'] = int(58)
dataByStim.loc[(dataByStim.index.str.find("A2") != -1)
               | (dataByStim.index.str.find("B2") != -1), 'AmbientLAeq'] = int(52)

# UAS SNR category
dataByStim.loc[(dataByStim.CALBINRecFiles.str.find("B2_CALBIN") != -1)
               | (dataByStim.CALBINRecFiles.str.find("A1_CALBIN") != -1)
               | (dataByStim.CALBINRecFiles.str.find("A2_CALBIN") != -1),
               'SNRlevel'] = "No UAS"
mask = (((dataByStim.index.str.find("A1") != -1)
         | (dataByStim.index.str.find("A2") != -1)
         | (dataByStim.index.str.find("B2") != -1))
        & ~((dataByStim.CALBINRecFiles.str.find("A1_CALBIN") != -1)
            | (dataByStim.CALBINRecFiles.str.find("A2_CALBIN") != -1)
            | (dataByStim.CALBINRecFiles.str.find("B2_CALBIN") != -1)))

dataByStim.loc[mask, 'SNRlevel'] = (dataByStim.loc[mask, 'UASLAeq']
                                    - dataByStim.loc[mask, 'AmbientLAeq']).astype(int)

# UAS operation
dataByStim.loc[dataByStim.index.str.find("_F_") != -1, 'UASOperation'] = "Flyover"
dataByStim.loc[dataByStim.index.str.find("_T_") != -1, 'UASOperation'] = "Takeoff"
dataByStim.loc[(dataByStim.index.str.find("_L_") != -1)
               & (dataByStim.index.str.find("_F_") == -1),
               'UASOperation'] = "Landing"
dataByStim.loc[(dataByStim.CALBINRecFiles.str.find("B2_CALBIN") != -1)
               | (dataByStim.CALBINRecFiles.str.find("A1_CALBIN") != -1)
               | (dataByStim.CALBINRecFiles.str.find("A2_CALBIN") != -1),
               'UASOperation'] = "No UAS"

# UAS event quantity / density
dataByStim.loc[dataByStim['UASOperation'] == "No UAS", 'UASEvents'] = 0
dataByStim.loc[(dataByStim['SessionPart'] == "A")
               & (dataByStim['UASOperation'] != "No UAS"), 'UASEvents'] = 1
dataByStim.loc[(dataByStim['SessionPart'] == "B")
               & (dataByStim.CALBINRecFiles.str.find("_1_CALBIN") != -1),
               'UASEvents'] = 1
dataByStim.loc[(dataByStim['SessionPart'] == "B")
               & (dataByStim.CALBINRecFiles.str.find("_3_CALBIN") != -1),
               'UASEvents'] = 3
dataByStim.loc[(dataByStim['SessionPart'] == "B")
               & (dataByStim.CALBINRecFiles.str.find("_5_CALBIN") != -1),
               'UASEvents'] = 5
dataByStim.loc[(dataByStim['SessionPart'] == "B")
               & (dataByStim.CALBINRecFiles.str.find("_9_CALBIN") != -1),
               'UASEvents'] = 9

dataByStim['UASEventPerMin'] = dataByStim['UASEvents']/(dataByStim['StimDuration']/60)


# UAS type
dataByStim.loc[dataByStim.index.str.find("_H520_") != -1, 'UASType'] = "H520"
dataByStim.loc[dataByStim.index.str.find("_M300_") != -1, 'UASType'] = "M300"
dataByStim.loc[dataByStim.index.str.find("_T150_") != -1, 'UASType'] = "T150"
dataByStim.loc[(dataByStim.CALBINRecFiles.str.find("B2_CALBIN") != -1)
               | (dataByStim.CALBINRecFiles.str.find("A1_CALBIN") != -1)
               | (dataByStim.CALBINRecFiles.str.find("A2_CALBIN") != -1),
               'UASType'] = "No UAS"

# ambient sound environment
dataByStim.loc[(dataByStim.index.str.find("A2") != -1)
               | (dataByStim.index.str.find("B2") != -1),
               'AmbientEnv'] = "Park"
dataByStim.loc[dataByStim.index.str.find("A1") != -1, 'AmbientEnv'] = "Street"

# -----------------------------------------------
# Acoustic, PNL & Detection metrics single values
# -----------------------------------------------


# output variables
indicesAcoustic = ["LAeqMaxLR", "LAEMaxLR", "LAFmaxMaxLR",
                   "LAF5ExMaxLR", "LAF10ExMaxLR", "LAF25ExMaxLR",
                   "LAF50ExMaxLR", "LAF75ExMaxLR", "LAF90ExMaxLR",
                   "LAF95ExMaxLR", "LASmaxMaxLR"]

dataByStim = pd.concat([dataByStim, pd.DataFrame(index=dataByStim.index,
                                                 columns=indicesAcoustic,
                                                 dtype=float)], axis=1)

# Acoustic metrics calculations
# -----------------------------

# loop over files to analyse
for ii, file in enumerate(filelist):
    if ii == 0:
        print("Processing acoustic metrics...\n")
    print(file.split('/')[-1])
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

    # take max L/R ear only
    dataByStim.loc[stemNames[ii], 'LAeqMaxLR'] = signalLAeq.max()
    dataByStim.loc[stemNames[ii], 'LAEMaxLR'] = signalLAE.max()
    dataByStim.loc[stemNames[ii], 'LAFmaxMaxLR'] = signalLAFmax.max()
    dataByStim.loc[stemNames[ii], 'LAF5ExMaxLR'] = signalLAF5.max()
    dataByStim.loc[stemNames[ii], 'LAF10ExMaxLR'] = signalLAF10.max()
    dataByStim.loc[stemNames[ii], 'LAF25ExMaxLR'] = signalLAF25.max()
    dataByStim.loc[stemNames[ii], 'LAF50ExMaxLR'] = signalLAF50.max()
    dataByStim.loc[stemNames[ii], 'LAF75ExMaxLR'] = signalLAF75.max()
    dataByStim.loc[stemNames[ii], 'LAF90ExMaxLR'] = signalLAF90.max()
    dataByStim.loc[stemNames[ii], 'LAF95ExMaxLR'] = signalLAF95.max()
    dataByStim.loc[stemNames[ii], 'LASmaxMaxLR'] = signalLASmax.max()

    # calculate intermittency ratio for test sounds
    C = [2, 3, 5]  # define threshold constant (dB)

    if (stemNames[ii].find("A1") != -1 or stemNames[ii].find("A2") != -1
        or stemNames[ii].find("B2") != -1):

        # 1-second LAeq
        signalLAeq1s = 20*np.log10(np.sqrt(pd.DataFrame(signalA[start_skips:-end_skips]**2).rolling(window=sampleRatein,
                                                                                                    step=sampleRatein,
                                                                                                    closed='left').mean())/2e-5)
        intermitRatioAll = np.zeros((len(C), 2))
        for jj, jjC in enumerate(C):
            intermitRatioK = signalLAeq + jjC  # IR threshold
            mask = signalLAeq1s > intermitRatioK
            # calculate average intensity for events exceeding threshold
            IAeqEvents = np.nansum(10**(signalLAeq1s[mask]/10), axis=0)/(len(signalA[start_skips:-end_skips])/sampleRatein)
            # convert to sound pressure while avoiding log(0)
            LAeqEvents = 10*np.log10(IAeqEvents, out=np.zeros_like(IAeqEvents, dtype=np.float64),
                                     where=(IAeqEvents != 0))
            intermitRatio = 100*10**((LAeqEvents - signalLAeq)/10)
            intermitRatio[IAeqEvents == 0] = 0  # replace tiny decimal with 0
            
            intermitRatioAll[jj, :] = intermitRatio
            
            dataByStim.loc[stemNames[ii], ('IntermitRatioC' + str(jjC) + 'MaxLR')] = intermitRatioAll[jj, :].max()




# end of for loop over CALBIN signal wav files


# PNL and detection metrics import
# --------------------------------

indicesPNL = ["PNLmaxMaxLR", "PNLTmaxMaxLR", "EPNLMaxLR"]
indicesDetect = ["Detect0p5dBADiscMaxLR", "Detect0p5dBMaxMaxLR",
                 "Detect0p5dBIntMaxLR", "Detect0p1dBADiscMaxLR",
                 "Detect0p1dBMaxMaxLR", "Detect0p1dBIntMaxLR"]

# import PNL and detection results
fileExts = "*.csv"
filelist = list(QFileDialog.getOpenFileNames(caption="Open results files in '03 Experiment\Experiment 1\Analysis\MATLAB\CALBIN\csv'",
                                             filter=fileExts))[0]
# assumes read files in alphabetical order with NASA_ first and SQAT_ second
DetectResults = pd.read_csv(filelist[0], header=0, index_col=0)
PNLResults = pd.read_csv(filelist[1], header=0, index_col=0)
PNLDetectResults = PNLResults.merge(DetectResults, left_index=True, right_index=True)
PNLDetectResults.set_index(PNLDetectResults.index.str.replace(pat="_CALBIN_Pa.wav", repl=""), inplace=True)

# join results to data
dataByStim = dataByStim.join(PNLDetectResults, how='left')


# ------------------------------------------------
# Psychoacoustic metrics single values calculation
# ------------------------------------------------

# output variables
indicesPsycho = ["LoudQZ5321PowAvgMaxLR",
                 "LoudQZ532105ExMaxLR",
                 "LoudQZ226PowAvgMaxLR",
                 "LoudQZ22605ExMaxLR",
                 "LoudQZ4182PowAvgMaxLR",
                 "LoudQZ418205ExMaxLR",
                 "LoudECMAPowAvgBin",
                 "SharpAurQZ5321PowAvgMaxLR",
                 "SharpAurQZ532105ExMaxLR",
                 "SharpAurQZ226PowAvgMaxLR",
                 "SharpAurQZ22605ExMaxLR",
                 "SharpAurQZ4182PowAvgMaxLR",
                 "SharpAurQZ418205ExMaxLR",
                 "LoudISO105ExMaxLR",
                 "LoudISO1PowAvgMaxLR",
                 "LoudISO3PowAvgBin",
                 "TonalECMAAvgMaxLR",
                 "TonalECMA05ExMaxLR",
                 "TonalSHMIntAvgMaxLR",
                 "TonalSHMInt05ExMaxLR",
                 "TonLdECMAPowAvgBin",
                 "TonLdECMA05ExBin",
                 "TonalAurAvgMaxLR",
                 "TonalAur10ExMaxLR",
                 "TonalAur05ExMaxLR",
                 "RoughECMA10ExBin",
                 "RoughECMA05ExBin",
                 "RoughFZ10ExMaxLR",
                 "RoughFZ05ExMaxLR",
                 "RoughDW10ExMaxLR",
                 "RoughDW05ExMaxLR",
                 "FluctSHM10ExBin",
                 "FluctSHM05ExBin",
                 "FluctFZ10ExMaxLR",
                 "FluctFZ05ExMaxLR",
                 "FluctOV10ExMaxLR",
                 "FluctOV05ExMaxLR",
                 "SharpAurSHMPowAvgMaxLR",
                 "SharpAurSHM05ExMaxLR",
                 "SharpAurISO3PowAvgMaxLR",
                 "SharpAurISO305ExMaxLR",
                 "SharpAurISO1PowAvgMaxLR",
                 "SharpAurISO105ExMaxLR",
                 "SharpDINPowAvgMaxLR",
                 "SharpDIN05ExMaxLR",
                 "SharpvBISO1PowAvgMaxLR",
                 "SharpvBISO105ExMaxLR",
                 "ImpulsSHMPowAvgMaxLR",
                 "ImpulsSHMAvgMaxLR",
                 "ImpulsSHM05ExMaxLR",
                 "dTonalECMAAvgMaxLR",
                 "dTonalECMA05ExMaxLR",
                 "dTonalSHMIntAvgMaxLR",
                 "dTonalSHMInt05ExMaxLR",
                 "dTonLdECMAPowAvgBin",
                 "dTonLdECMA05ExBin",
                 "dRoughECMA10ExBin",
                 "dRoughECMA05ExBin",
                 "dRoughFZ10ExMaxLR",
                 "dRoughFZ05ExMaxLR",
                 "dFluctSHM10ExBin",
                 "dFluctSHM05ExBin",
                 "dFluctOV10ExMaxLR",
                 "dFluctOV05ExMaxLR",
                 "dSharpAurSHMPowAvgMaxLR",
                 "dSharpAurSHM05ExMaxLR",
                 "dSharpAurISO3PowAvgMaxLR",
                 "dSharpAurISO305ExMaxLR",
                 "dImpulsSHMPowAvgMaxLR",
                 "dImpulsSHMAvgMaxLR",
                 "dImpulsSHM05ExMaxLR"]

dataByStim = pd.concat([dataByStim, pd.DataFrame(index=dataByStim.index,
                                                 columns=indicesPsycho,
                                                 dtype=float)], axis=1)

# output SQM sample rates
sampleRateLoudECMA = 187.5
sampleRateLoudISO1 = 500
sampleRateLoudISO3 = 1e3
sampleRateTonalECMA = 187.5
sampleRateTonalAures = 12.5
sampleRateRoughECMA = 50
sampleRateRoughFZ = 2000
sampleRateRoughDW = 10
sampleRateFluctSHM = 229390681/4e6
sampleRateFluctOV = 5
sampleRateImpulsSHM = 60000/55
sampleRateQZ = 10

# critical band differences (used for integration of specific values)
bandDiff0p5 = 0.5  # used for all except Sottek Hearing Model impulsiveness
bandDiffSHMImp = 1.0  # Sottek Hearing Model impulsiveness band differences

# time value for moving averaging of SQMs for difference calculations
windowT = 0.05

# From MATLAB calculations
# ------------------------

# open xlsx file selection dialog and assign filepaths to list
# PROJECT NOTE: the calculated files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Analysis/MATLAB/CALHEQ?csf=1&web=1&e=BNhgIT
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

fileExtsXLSX = "*.xlsx"
filelistXLSX = list(QFileDialog.getOpenFileNames(filter=fileExtsXLSX,
                                                 caption=r"Select MATLAB output files in '03 Experiment\Experiment 1\Analysis\MATLAB\CALHEQ\xlsx'"))[0]
filelistXLSX.sort()

fileExtsMAT = "*.mat"
filelistMAT = list(QFileDialog.getOpenFileNames(filter=fileExtsMAT,
                                                caption=r"Select MATLAB output files in '03 Experiment\Experiment 1\Analysis\MATLAB\CALHEQ\mat'"))[0]

filelistMAT.sort()

filelist = filelistXLSX + filelistMAT

filenames = [filepath.split('/')[-1] for filepath in filelist]
renderNames = [filename.split('_CALHEQ_Pa')[0] for filename in filenames]

# loop over files to analyse
for ii, file in enumerate(filelist):
    if ii == 0:
        print("Processing sound quality metrics...\n")
    print(file.split('/')[-1])

    if file[-3:] == "mat":
        if filenames[ii].find("Fluct") != -1:
            try:
                filedata = io.loadmat(file)['FluctFZSpecTDep2']
            except:
                filedata = io.loadmat(file)['FluctFZSpecTDep']

            chan2col = np.arange(0, filedata.shape[1])[filedata[0, :] == 0.5][-1]

            # Calculate overall Fastl & Zwicker fluctuation strength from specific
            # time-dependent fluctuation strength
            specFluctFZL = pd.DataFrame(filedata[1:, 1:chan2col],
                                        index=filedata[1:, 0],
                                        columns=filedata[0, 1:chan2col])
            specFluctFZR = pd.DataFrame(filedata[1:, chan2col:],
                                        index=filedata[1:, 0],
                                        columns=filedata[0, chan2col:])

            fluctFZTDep = pd.merge(left=bandDiff0p5*specFluctFZL.sum(axis=1).to_frame(name="Channel 1"),
                                   right=bandDiff0p5*specFluctFZR.sum(axis=1).to_frame(name="Channel 2"),
                                   left_index=True, right_index=True)

            # mask for start/end skip and 0 values
            fluctFZTDepMask = fluctFZTDep.loc[(fluctFZTDep.index.values
                                               > start_skipT).transpose()
                                              & (fluctFZTDep.index.values
                                                 < fluctFZTDep.index.values.max()
                                                 - end_skipT).transpose()]

            # 2-channel overall 5% exceeded fluctuation strength
            fluctFZ05Ex = fluctFZTDepMask.quantile(q=0.95)
            # max of l/r channel overall 5% exceeded fluctuation strength
            fluctFZ05ExMaxLR = fluctFZ05Ex.max()

            # 2-channel overall 10% exceeded fluctuation strength
            fluctFZ10Ex = fluctFZTDepMask.quantile(q=0.90)
            # max of l/r channel overall 10% exceeded fluctuation strength
            fluctFZ10ExMaxLR = fluctFZ10Ex.max()
            
            dataByStim.loc[renderNames[ii], 'FluctFZ10ExMaxLR'] = fluctFZ10ExMaxLR
            dataByStim.loc[renderNames[ii], 'FluctFZ05ExMaxLR'] = fluctFZ05ExMaxLR

        elif filenames[ii].find("Rough") != -1:
            try:
                filedata = io.loadmat(file)['RoughFZSpecTDep2']
            except:
                filedata = io.loadmat(file)['RoughFZSpecTDep']

            chan2col = np.arange(0, filedata.shape[1])[filedata[0, :] == 0.5][-1]

            # Calculate overall Fastl & Zwicker roughness from specific
            # time-dependent roughness
            specRoughFZL = pd.DataFrame(filedata[1:, 1:chan2col],
                                        index=filedata[1:, 0],
                                        columns=filedata[0, 1:chan2col])
            specRoughFZR = pd.DataFrame(filedata[1:, chan2col:],
                                        index=filedata[1:, 0],
                                        columns=filedata[0, chan2col:])

            roughFZTDep = pd.merge(left=bandDiff0p5*specRoughFZL.sum(axis=1).to_frame(name="Channel 1"),
                                   right=bandDiff0p5*specRoughFZR.sum(axis=1).to_frame(name="Channel 2"),
                                   left_index=True, right_index=True)

            # mask for start/end skip and 0 values
            roughFZTDepMask = roughFZTDep.loc[(roughFZTDep.index.values
                                               > start_skipT).transpose()
                                              & (roughFZTDep.index.values
                                                 < roughFZTDep.index.values.max()
                                                 - end_skipT).transpose()]

            # 2-channel overall 5% exceeded roughness
            roughFZ05Ex = roughFZTDepMask.quantile(q=0.95)
            # max of l/r channel overall 5% exceeded roughness
            roughFZ05ExMaxLR = roughFZ05Ex.max()

            # 2-channel overall 10% exceeded roughness
            roughFZ10Ex = roughFZTDepMask.quantile(q=0.90)
            # max of l/r channel overall 10% exceeded roughness
            roughFZ10ExMaxLR = roughFZ10Ex.max()

            # add results to output DataFrame
            # dataByStim.loc[renderNames[ii], 'RoughFZ10ExMaxLR'] = roughFZ10ExMaxLR
            # dataByStim.loc[renderNames[ii], 'RoughFZ05ExMaxLR'] = roughFZ05ExMaxLR

            # calculation section for SQM differences
            # NOTE: THIS SECTION RELIES ON THE ALPHABETIC ORDER OF THE STIMULI FILES AS
            # ORIGINALLY NAMED: EACH AMBIENT SOUND FILE PRECEDING THE CORRESPONDING
            # COMBINED STIMULI FILES. IF THE ORDERING OR FILE NAMING IS CHANGED, THIS
            # CALCULATION WILL BE INVALIDATED AND *MAY ALSO* CAUSE AN ERROR.
            # a rolling 50 ms window is applied to average the SQM values over time -
            # this is to reduce uncertainty due to imperfect time-alignment between the
            # ambient vs combined stimuli files (all are from recordings, so there will
            # be some slippage due to imperfect editing)
            if renderNames[ii] in ["A1", "A2", "B2"]:
                # NOTE: we could dropna() the first <windowT values, but these will be
                # ignored anyway in the statistical analysis, assuming start_skipT >
                # windowT

                # calculate moving average values for ambient stimulus
                ambSpecRoughFZLMovAvg = specRoughFZL.rolling(window=int(np.ceil(sampleRateRoughFZ*windowT))).mean().T
                ambSpecRoughFZRMovAvg = specRoughFZR.rolling(window=int(np.ceil(sampleRateRoughFZ*windowT))).mean().T

            elif renderNames[ii][0:3] in ["A1_", "A2_", "B2_"]:
                # calculate moving average values for combined stimulus
                specRoughFZLMovAvg = specRoughFZL.rolling(window=int(np.ceil(sampleRateRoughFZ*windowT))).mean().T
                specRoughFZRMovAvg = specRoughFZR.rolling(window=int(np.ceil(sampleRateRoughFZ*windowT))).mean().T

                # calculate differences and make negative values 0
                dSpecRoughFZL = np.maximum(specRoughFZLMovAvg
                                           - ambSpecRoughFZLMovAvg, 0)
                dSpecRoughFZR = np.maximum(specRoughFZRMovAvg
                                           - ambSpecRoughFZRMovAvg, 0)

                # calculate aggregated difference values

                # 2-channel time-dependent roughness
                dRoughFZTDep = pd.concat([bandDiff0p5*dSpecRoughFZL.sum(axis=0),
                                          bandDiff0p5*dSpecRoughFZR.sum(axis=0)],
                                         axis=1)

                # mask for start/end skip
                dRoughFZTDepMask = roughFZTDep.loc[(dRoughFZTDep.index.values
                                                    > start_skipT)
                                                   & (dRoughFZTDep.index.values
                                                      < dRoughFZTDep.index.values.max()
                                                      - end_skipT)]

                # overall (90th percentile = 10% exceeded) roughness
                dRoughFZ10Ex = dRoughFZTDepMask.quantile(q=0.90)
                # max of l/r channel overall 10% exceeded roughness
                dRoughFZ10ExMaxLR = dRoughFZ10Ex.max()
                # overall (95th percentile = 5% exceeded) roughness
                dRoughFZ05Ex = dRoughFZTDepMask.quantile(q=0.95)
                # max of l/r channel overall 5% exceeded roughness
                dRoughFZ05ExMaxLR = dRoughFZ05Ex.max()

                # add results to output DataFrame
                dataByStim.loc[renderNames[ii], 'dRoughFZ10ExMaxLR'] = dRoughFZ10ExMaxLR
                dataByStim.loc[renderNames[ii], 'dRoughFZ05ExMaxLR'] = dRoughFZ05ExMaxLR

    # end of if branch for .mat files

    elif file[-4:] == "xlsx":
        workbookdata = pd.read_excel(io=file, sheet_name=None)

        # Calculate quasi-Zwicker (ISO 532-1) overall loudness from 2-channel
        # time-dependent loudness
        loudQZ5321TDep = pd.DataFrame(workbookdata['LoudQZ5321'].iloc[0:, 1:3].values,
                                         columns=workbookdata['LoudQZ5321'].iloc[0, 1:3].index,
                                         index=workbookdata['LoudQZ5321'].iloc[:, 0])

        # mask for start/end skip
        loudQZ5321TDepMask = loudQZ5321TDep.loc[(loudQZ5321TDep.index.values
                                               > start_skipT)
                                              & (loudQZ5321TDep.index.values
                                                 < loudQZ5321TDep.index.values.max()
                                                 - end_skipT), :]

        # 2-channel overall (power-averaged) loudness
        loudQZ5321PowAvg = ((loudQZ5321TDepMask**(1/np.log10(2))).sum(axis=0)
                           / len(loudQZ5321TDepMask))**np.log10(2)
        # max of l/r channel overall (power-averaged) loudness
        loudQZ5321PowAvgMaxLR = loudQZ5321PowAvg.max()
        # 2-channel overall 5% exceeded loudness
        loudQZ532105Ex = loudQZ5321TDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded loudness
        loudQZ532105ExMaxLR = loudQZ532105Ex.max()
        
        # Calculate overall Aures+quasi-Zwicker (ISO 226 adjusted) loudness from 2-channel
        # time-dependent loudness
        loudQZ226TDep = pd.DataFrame(workbookdata['LoudQZ226'].iloc[0:, 1:3].values,
                                         columns=workbookdata['LoudQZ226'].iloc[0, 1:3].index,
                                         index=workbookdata['LoudQZ226'].iloc[:, 0])

        # mask for start/end skip
        loudQZ226TDepMask = loudQZ226TDep.loc[(loudQZ226TDep.index.values
                                               > start_skipT)
                                              & (loudQZ226TDep.index.values
                                                 < loudQZ226TDep.index.values.max()
                                                 - end_skipT), :]

        # 2-channel overall (power-averaged) loudness
        loudQZ226PowAvg = ((loudQZ226TDepMask**(1/np.log10(2))).sum(axis=0)
                           / len(loudQZ226TDepMask))**np.log10(2)
        # max of l/r channel overall (power-averaged) loudness
        loudQZ226PowAvgMaxLR = loudQZ226PowAvg.max()
        # 2-channel overall 5% exceeded loudness
        loudQZ22605Ex = loudQZ226TDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded loudness
        loudQZ22605ExMaxLR = loudQZ22605Ex.max()
        
        # Calculate overall Aures+quasi-Zwicker (ECMA-418-2 ear filtered) loudness from 2-channel
        # time-dependent loudness
        loudQZ4182TDep = pd.DataFrame(workbookdata['LoudQZ4182'].iloc[0:, 1:3].values,
                                         columns=workbookdata['LoudQZ4182'].iloc[0, 1:3].index,
                                         index=workbookdata['LoudQZ4182'].iloc[:, 0])

        # mask for start/end skip
        loudQZ4182TDepMask = loudQZ4182TDep.loc[(loudQZ4182TDep.index.values
                                               > start_skipT)
                                              & (loudQZ4182TDep.index.values
                                                 < loudQZ4182TDep.index.values.max()
                                                 - end_skipT), :]

        # 2-channel overall (power-averaged) loudness
        loudQZ4182PowAvg = ((loudQZ4182TDepMask**(1/np.log10(2))).sum(axis=0)
                           / len(loudQZ4182TDepMask))**np.log10(2)
        # max of l/r channel overall (power-averaged) loudness
        loudQZ4182PowAvgMaxLR = loudQZ4182PowAvg.max()
        # 2-channel overall 5% exceeded loudness
        loudQZ418205Ex = loudQZ4182TDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded loudness
        loudQZ418205ExMaxLR = loudQZ418205Ex.max()
        
        # Calculate overall Aures+quasi-Zwicker sharpness from 2-channel
        # time-dependent sharpness
        sharpAQZ5321TDep = pd.DataFrame(workbookdata['SharpAuresQZ5321'].iloc[0:, 1:3].values,
                                         columns=workbookdata['SharpAuresQZ5321'].iloc[0, 1:3].index,
                                         index=workbookdata['SharpAuresQZ5321'].iloc[:, 0])

        # mask for start/end skip
        sharpAQZ5321TDepMask = sharpAQZ5321TDep.loc[(sharpAQZ5321TDep.index.values
                                               > start_skipT)
                                              & (sharpAQZ5321TDep.index.values
                                                 < sharpAQZ5321TDep.index.values.max()
                                                 - end_skipT), :]

        # 2-channel overall (power-averaged) sharpness
        sharpAQZ5321PowAvg = ((sharpAQZ5321TDepMask**(1/np.log10(2))).sum(axis=0)
                           / len(sharpAQZ5321TDepMask))**np.log10(2)
        # max of l/r channel overall (power-averaged) sharpness
        sharpAQZ5321PowAvgMaxLR = sharpAQZ5321PowAvg.max()
        # 2-channel overall 5% exceeded sharpness
        sharpAQZ532105Ex = sharpAQZ5321TDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded sharpness
        sharpAQZ532105ExMaxLR = sharpAQZ532105Ex.max()
        
        # Calculate overall Aures+quasi-Zwicker (ISO 226 adjusted) sharpness from 2-channel
        # time-dependent sharpness
        sharpAQZ226TDep = pd.DataFrame(workbookdata['SharpAuresQZ226'].iloc[0:, 1:3].values,
                                         columns=workbookdata['SharpAuresQZ226'].iloc[0, 1:3].index,
                                         index=workbookdata['SharpAuresQZ226'].iloc[:, 0])

        # mask for start/end skip
        sharpAQZ226TDepMask = sharpAQZ226TDep.loc[(sharpAQZ226TDep.index.values
                                               > start_skipT)
                                              & (sharpAQZ226TDep.index.values
                                                 < sharpAQZ226TDep.index.values.max()
                                                 - end_skipT), :]

        # 2-channel overall (power-averaged) sharpness
        sharpAQZ226PowAvg = ((sharpAQZ226TDepMask**(1/np.log10(2))).sum(axis=0)
                           / len(sharpAQZ226TDepMask))**np.log10(2)
        # max of l/r channel overall (power-averaged) sharpness
        sharpAQZ226PowAvgMaxLR = sharpAQZ226PowAvg.max()
        # 2-channel overall 5% exceeded sharpness
        sharpAQZ22605Ex = sharpAQZ226TDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded sharpness
        sharpAQZ22605ExMaxLR = sharpAQZ22605Ex.max()
        
        # Calculate overall Aures+quasi-Zwicker (ECMA-418-2 ear filtered) sharpness from 2-channel
        # time-dependent sharpness
        sharpAQZ4182TDep = pd.DataFrame(workbookdata['SharpAuresQZ4182'].iloc[0:, 1:3].values,
                                         columns=workbookdata['SharpAuresQZ4182'].iloc[0, 1:3].index,
                                         index=workbookdata['SharpAuresQZ4182'].iloc[:, 0])

        # mask for start/end skip
        sharpAQZ4182TDepMask = sharpAQZ4182TDep.loc[(sharpAQZ4182TDep.index.values
                                               > start_skipT)
                                              & (sharpAQZ4182TDep.index.values
                                                 < sharpAQZ4182TDep.index.values.max()
                                                 - end_skipT), :]

        # 2-channel overall (power-averaged) sharpness
        sharpAQZ4182PowAvg = ((sharpAQZ4182TDepMask**(1/np.log10(2))).sum(axis=0)
                              / len(sharpAQZ4182TDepMask))**np.log10(2)
        # max of l/r channel overall (power-averaged) sharpness
        sharpAQZ4182PowAvgMaxLR = sharpAQZ4182PowAvg.max()
        # 2-channel overall 5% exceeded sharpness
        sharpAQZ418205Ex = sharpAQZ4182TDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded sharpness
        sharpAQZ418205ExMaxLR = sharpAQZ418205Ex.max()

        # Calculate ECMA-418-2:2022 Sottek Hearing Model overall loudness from
        # 2-channel specific loudness
        # left channel
        specLoudECMAL = pd.DataFrame(workbookdata['LoudECMAL'].iloc[:, 1:].values.T,
                                      columns=workbookdata['LoudECMAL'].iloc[:, 0],
                                      index=workbookdata['LoudECMAL'].iloc[0, 1:].index)
        # right channel
        specLoudECMAR = pd.DataFrame(workbookdata['LoudECMAR'].iloc[:, 1:].values.T,
                                      columns=workbookdata['LoudECMAR'].iloc[:, 0],
                                      index=workbookdata['LoudECMAR'].iloc[0, 1:].index)

        # binaural specific loudness (ECMA-418-2:2022 Equation 118)
        specLoudECMABin = ((specLoudECMAL**2
                            + specLoudECMAR**2)/2).pow(0.5)

        # binaural time-dependent loudness (ECMA-418-2:2022 Equation 116)
        loudECMATDepBin = specLoudECMABin.sum(axis=0)*bandDiff0p5

        # mask for start/end skip
        loudECMATDepBinMask = loudECMATDepBin.loc[(loudECMATDepBin.index.values
                                                    > start_skipT)
                                                  & (loudECMATDepBin.index.values
                                                      < loudECMATDepBin.index.values.max()
                                                      - end_skipT)]

        # binaural overall (power-averaged) loudness (ECMA-418-2:2022 Equation 117)
        loudECMAPowAvgBin = ((loudECMATDepBinMask**(1/np.log10(2))).sum()
                              / len(loudECMATDepBinMask))**np.log10(2)

        # Calculate ECMA-418-2:2022 Sottek Hearing Model overall tonality from
        # 2-channel specific tonality
        # left channel
        specTonalECMAL = pd.DataFrame(workbookdata['TonalECMAL'].iloc[:, 1:].values.T,
                                      columns=workbookdata['TonalECMAL'].iloc[:, 0],
                                      index=workbookdata['TonalECMAL'].iloc[0, 1:].index)
        # right channel
        specTonalECMAR = pd.DataFrame(workbookdata['TonalECMAR'].iloc[:, 1:].values.T,
                                      columns=workbookdata['TonalECMAR'].iloc[:, 0],
                                      index=workbookdata['TonalECMAR'].iloc[0, 1:].index)

        # 2-channel time-dependent tonality (max, not integration)
        tonalECMATDep = pd.concat([specTonalECMAL.max(axis=0),
                                    specTonalECMAR.max(axis=0)],
                                  axis=1)

        # 2-channel time-dependent tonality (integrated with adjustment to match 40 dB 1 kHz sine to 1 tu)
        tonalSHMIntTDep = pd.concat([specTonalECMAL.sum(axis=0)*bandDiff0p5*0.348088948583815,
                                      specTonalECMAR.sum(axis=0)*bandDiff0p5*0.348088948583815],
                                    axis=1)

        # mask for start/end skip and values <= 0.02
        tonalECMATDepMaskL = tonalECMATDep.loc[(tonalECMATDep.index.values
                                                > start_skipT)
                                                & (tonalECMATDep.index.values
                                                  < tonalECMATDep.index.values.max()
                                                  - end_skipT)
                                                & (tonalECMATDep.loc[:, 0].values
                                                  > 0.02), 0]
        tonalECMATDepMaskR = tonalECMATDep.loc[(tonalECMATDep.index.values
                                                > start_skipT)
                                                & (tonalECMATDep.index.values
                                                  < tonalECMATDep.index.values.max()
                                                  - end_skipT)
                                                & (tonalECMATDep.loc[:, 1].values
                                                  > 0.02), 1]

        # mask for start/end skip and values <= 0.02
        # NOTE: uses mask from ECMA tonality
        tonalSHMIntTDepMaskL = tonalSHMIntTDep.loc[(tonalSHMIntTDep.index.values
                                                    > start_skipT)
                                                    & (tonalSHMIntTDep.index.values
                                                      < tonalSHMIntTDep.index.values.max()
                                                      - end_skipT)
                                                    & (tonalECMATDep.loc[:, 0].values
                                                      > 0.02), 0]  # see NOTE above
        tonalSHMIntTDepMaskR = tonalSHMIntTDep.loc[(tonalSHMIntTDep.index.values
                                                    > start_skipT)
                                                    & (tonalSHMIntTDep.index.values
                                                      < tonalSHMIntTDep.index.values.max()
                                                      - end_skipT)
                                                    & (tonalECMATDep.loc[:, 1].values
                                                      > 0.02), 1]  # see NOTE above

        # 2-channel time-averaged tonality (omitting T<=0.02)
        tonalECMAAvgL = tonalECMATDepMaskL.mean(axis=0)
        tonalECMAAvgR = tonalECMATDepMaskR.mean(axis=0)
        # max of L/R
        tonalECMAAvgMaxLR = max(tonalECMAAvgL, tonalECMAAvgR)

        # 2-channel 5% exceeded tonality (omitting T<=0.02)
        tonalECMA05ExL = tonalECMATDepMaskL.quantile(q=0.95)
        tonalECMA05ExR = tonalECMATDepMaskR.quantile(q=0.95)
        # max of L/R
        tonalECMA05ExMaxLR = max(tonalECMA05ExL, tonalECMA05ExR)

        # 2-channel time-averaged integated tonality (omitting T<=0.02)
        # NOTE: uses mask from ECMA tonality
        tonalSHMIntAvgL = tonalSHMIntTDepMaskL.mean(axis=0)
        tonalSHMIntAvgR = tonalSHMIntTDepMaskR.mean(axis=0)
        # max of L/R
        tonalSHMIntAvgMaxLR = max(tonalSHMIntAvgL, tonalSHMIntAvgR)

        # 2-channel 5% exceeded integrated tonality (omitting T<=0.02)
        tonalSHMInt05ExL = tonalSHMIntTDepMaskL.quantile(q=0.95)
        tonalSHMInt05ExR = tonalSHMIntTDepMaskR.quantile(q=0.95)
        # max of L/R
        tonalSHMInt05ExMaxLR = max(tonalSHMInt05ExL, tonalSHMInt05ExR)

        # Calculate ECMA-418-2:2022 Sottek Hearing Model overall tonal loudness from
        # 2-channel specific tonal loudness
        # left channel
        specTonLdECMAL = pd.DataFrame(workbookdata['TonLdECMAL'].iloc[:, 1:].values.T,
                                      columns=workbookdata['TonLdECMAL'].iloc[:, 0],
                                      index=workbookdata['TonLdECMAL'].iloc[0, 1:].index)
        # right channel
        specTonLdECMAR = pd.DataFrame(workbookdata['TonLdECMAR'].iloc[:, 1:].values.T,
                                      columns=workbookdata['TonLdECMAR'].iloc[:, 0],
                                      index=workbookdata['TonLdECMAR'].iloc[0, 1:].index)

        # binaural specific tonal loudness (ECMA-418-2:2022 Equation 118)
        specTonLdECMABin = ((specTonLdECMAL**2
                            + specTonLdECMAR**2)/2).pow(0.5)

        # binaural time-dependent tonal loudness (ECMA-418-2:2022 Equation 116)
        tonLdECMATDepBin = specTonLdECMABin.sum(axis=0)*bandDiff0p5

        # mask for start/end skip
        tonLdECMATDepBinMask = tonLdECMATDepBin.loc[(tonLdECMATDepBin.index.values
                                                      > start_skipT)
                                                    & (tonLdECMATDepBin.index.values
                                                        < tonLdECMATDepBin.index.values.max()
                                                        - end_skipT)]

        # binaural overall (power-averaged) tonal loudness (ECMA-418-2:2022
        # Equation 117)
        tonLdECMAPowAvgBin = ((tonLdECMATDepBinMask**(1/np.log10(2))).sum()
                              / len(tonLdECMATDepBinMask))**np.log10(2)

        # binaural 5% exceeded tonal loudness
        tonLdECMA05ExBin = tonLdECMATDepBinMask.quantile(q=0.95)

        # Calculate overall Aures tonality from 2-channel time-dependent
        # tonality
        tonalAurTDep = pd.DataFrame(workbookdata['TonalAures'].iloc[0:, 1:3].values,
                                    columns=workbookdata['TonalAures'].iloc[0, 1:3].index,
                                    index=workbookdata['TonalAures'].iloc[:, 0])

        # mask for start/end skip
        tonalAurTDepMask = tonalAurTDep.loc[(tonalAurTDep.index.values
                                              > start_skipT)
                                            & (tonalAurTDep.index.values
                                                < tonalAurTDep.index.values.max()
                                                - end_skipT), :]

        # 2-channel overall mean tonality
        tonalAurAvg = tonalAurTDepMask.mean(axis=0)
        # max of l/r channel overall mean tonality
        tonalAurAvgMaxLR = tonalAurAvg.max()
        # 2-channel overall 5% exceeded tonality
        tonalAur05Ex = tonalAurTDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded tonality
        tonalAur05ExMaxLR = tonalAur05Ex.max()
        # 2-channel overall 10% exceeded tonality
        tonalAur10Ex = tonalAurTDepMask.quantile(q=0.90)
        # max of l/r channel overall 10% exceeded tonality
        tonalAur10ExMaxLR = tonalAur10Ex.max()

        # Calculate ECMA-418-2:2022 Sottek Hearing Model binaural overall roughness
        # from specific roughness
        specRoughECMAL = pd.DataFrame(workbookdata['RoughECMAL'].iloc[:, 1:].values.T,
                                      columns=workbookdata['RoughECMAL'].iloc[:, 0],
                                      index=workbookdata['RoughECMAL'].iloc[0, 1:].index)

        specRoughECMAR = pd.DataFrame(workbookdata['RoughECMAR'].iloc[:, 1:].values.T,
                                      columns=workbookdata['RoughECMAR'].iloc[:, 0],
                                      index=workbookdata['RoughECMAR'].iloc[0, 1:].index)

        # binaural specific roughness (ECMA-418-2:2022 Equation 112)
        specRoughECMABin = ((specRoughECMAL**2
                            + specRoughECMAL**2)/2).pow(0.5)

        # binaural time-dependent roughness
        roughECMATDepBin = specRoughECMABin.sum(axis=0)*bandDiff0p5

        # mask for start/end skip
        roughECMATDepBinMask = roughECMATDepBin.loc[(roughECMATDepBin.index.values
                                                      > start_skipT)
                                                    & (roughECMATDepBin.index.values
                                                        < roughECMATDepBin.index.values.max()
                                                        - end_skipT)]

        # binaural overall (90th percentile = 10% exceeded) roughness
        roughECMA10ExBin = roughECMATDepBinMask.quantile(q=0.90)
        # binaural overall (95th percentile = 5% exceeded) roughness
        roughECMA05ExBin = roughECMATDepBinMask.quantile(q=0.95)

        # Calculate Daniel & Weber overall roughness from specific roughness
        specRoughDWL = pd.DataFrame(workbookdata['RoughDanWebL'].iloc[:, 1:].values.T,
                                    columns=workbookdata['RoughDanWebL'].iloc[:, 0],
                                    index=workbookdata['RoughDanWebL'].iloc[0, 1:].index)

        specRoughDWR = pd.DataFrame(workbookdata['RoughDanWebR'].iloc[:, 1:].values.T,
                                    columns=workbookdata['RoughDanWebR'].iloc[:, 0],
                                    index=workbookdata['RoughDanWebR'].iloc[0, 1:].index)

        # 2-channel time-dependent roughness
        roughDWTDep = pd.concat([bandDiff0p5*specRoughDWL.sum(axis=0),
                                  bandDiff0p5*specRoughDWR.sum(axis=0)],
                                axis=1)

        # mask for start/end skip
        roughDWTDepMask = roughDWTDep.loc[(roughDWTDep.index.values
                                            > start_skipT)
                                          & (roughDWTDep.index.values
                                              < roughDWTDep.index.values.max()
                                              - end_skipT)]

        # overall (90th percentile = 10% exceeded) roughness
        roughDW10Ex = roughDWTDepMask.quantile(q=0.90)
        # max of l/r channel overall 10% exceeded roughness
        roughDW10ExMaxLR = roughDW10Ex.max()
        # overall (95th percentile = 5% exceeded) roughness
        roughDW05Ex = roughDWTDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded roughness
        roughDW05ExMaxLR = roughDW05Ex.max()

        # Calculate Osses Vecchi et al overall fluctuation strength from
        # specific fluctuation strength
        specFluctOVL = pd.DataFrame(workbookdata['FluctOssVecL'].iloc[:, 1:].values.T,
                                    columns=workbookdata['FluctOssVecL'].iloc[:, 0],
                                    index=workbookdata['FluctOssVecL'].iloc[0, 1:].index)

        specFluctOVR = pd.DataFrame(workbookdata['FluctOssVecR'].iloc[:, 1:].values.T,
                                    columns=workbookdata['FluctOssVecR'].iloc[:, 0],
                                    index=workbookdata['FluctOssVecR'].iloc[0, 1:].index)

        # 2-channel time-dependent fluctuation strength
        fluctOVTDep = pd.concat([specFluctOVL.sum(axis=0)*bandDiff0p5,
                                  specFluctOVR.sum(axis=0)*bandDiff0p5],
                                axis=1)

        # mask for start/end skip
        fluctOVTDepMask = fluctOVTDep.loc[(fluctOVTDep.index.values
                                            > start_skipT)
                                          & (fluctOVTDep.index.values
                                              < fluctOVTDep.index.values.max()
                                              - end_skipT)]

        # overall (90th percentile = 10% exceeded) fluctuation strength
        fluctOV10Ex = fluctOVTDepMask.quantile(q=0.90)
        # max of l/r channel overall 10% exceeded fluctuation strength
        fluctOV10ExMaxLR = fluctOV10Ex.max()
        # overall (95th percentile = 5% exceeded) fluctuation strength
        fluctOV05Ex = fluctOVTDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded fluctuation strength
        fluctOV05ExMaxLR = fluctOV05Ex.max()

        # Calculate overall Aures+Sottek Hearing Model sharpness from 2-channel
        # time-dependent sharpness
        sharpASHMTDep = pd.DataFrame(workbookdata['SharpAuresSHM'].iloc[0:, 1:3].values,
                                      columns=workbookdata['SharpAuresSHM'].iloc[0, 1:3].index,
                                      index=workbookdata['SharpAuresSHM'].iloc[:, 0])

        # mask for start/end skip
        sharpASHMTDepMask = sharpASHMTDep.loc[(sharpASHMTDep.index.values
                                                > start_skipT)
                                              & (sharpASHMTDep.index.values
                                                  < sharpASHMTDep.index.values.max()
                                                  - end_skipT), :]

        # 2-channel overall (power-averaged) sharpness
        sharpASHMPowAvg = ((sharpASHMTDepMask**(1/np.log10(2))).sum(axis=0)
                            / len(sharpASHMTDepMask))**np.log10(2)
        # max of l/r channel overall (power-averaged) sharpness
        sharpASHMPowAvgMaxLR = sharpASHMPowAvg.max()
        # 2-channel overall 5% exceeded sharpness
        sharpASHM05Ex = sharpASHMTDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded sharpness
        sharpASHM05ExMaxLR = sharpASHM05Ex.max()

        # add results to output DataFrame
        dataByStim.loc[renderNames[ii], 'LoudQZ5321PowAvgMaxLR'] = loudQZ5321PowAvgMaxLR
        dataByStim.loc[renderNames[ii], 'LoudQZ532105ExMaxLR'] = loudQZ532105ExMaxLR
        dataByStim.loc[renderNames[ii], 'LoudQZ226PowAvgMaxLR'] = loudQZ226PowAvgMaxLR
        dataByStim.loc[renderNames[ii], 'LoudQZ22605ExMaxLR'] = loudQZ22605ExMaxLR
        dataByStim.loc[renderNames[ii], 'LoudQZ4182PowAvgMaxLR'] = loudQZ4182PowAvgMaxLR
        dataByStim.loc[renderNames[ii], 'LoudQZ418205ExMaxLR'] = loudQZ418205ExMaxLR
        dataByStim.loc[renderNames[ii], 'SharpAurQZ5321PowAvgMaxLR'] = sharpAQZ5321PowAvgMaxLR
        dataByStim.loc[renderNames[ii], 'SharpAurQZ532105ExMaxLR'] = sharpAQZ532105ExMaxLR
        dataByStim.loc[renderNames[ii], 'SharpAurQZ226PowAvgMaxLR'] = sharpAQZ226PowAvgMaxLR
        dataByStim.loc[renderNames[ii], 'SharpAurQZ22605ExMaxLR'] = sharpAQZ22605ExMaxLR
        dataByStim.loc[renderNames[ii], 'SharpAurQZ4182PowAvgMaxLR'] = sharpAQZ4182PowAvgMaxLR
        dataByStim.loc[renderNames[ii], 'SharpAurQZ418205ExMaxLR'] = sharpAQZ418205ExMaxLR
        dataByStim.loc[renderNames[ii], 'LoudECMAPowAvgBin'] = loudECMAPowAvgBin
        dataByStim.loc[renderNames[ii], 'TonalECMAAvgMaxLR'] = tonalECMAAvgMaxLR
        dataByStim.loc[renderNames[ii], 'TonalECMA05ExMaxLR'] = tonalECMA05ExMaxLR
        dataByStim.loc[renderNames[ii], 'TonalSHMIntAvgMaxLR'] = tonalSHMIntAvgMaxLR
        dataByStim.loc[renderNames[ii], 'TonalSHMInt05ExMaxLR'] = tonalSHMInt05ExMaxLR
        dataByStim.loc[renderNames[ii], 'TonLdECMAPowAvgBin'] = tonLdECMAPowAvgBin
        dataByStim.loc[renderNames[ii], 'TonLdECMA05ExBin'] = tonLdECMA05ExBin
        dataByStim.loc[renderNames[ii], 'TonalAur05ExMaxLR'] = tonalAur05ExMaxLR
        dataByStim.loc[renderNames[ii], 'TonalAur10ExMaxLR'] = tonalAur10ExMaxLR
        dataByStim.loc[renderNames[ii], 'TonalAurAvgMaxLR'] = tonalAurAvgMaxLR
        dataByStim.loc[renderNames[ii], 'RoughECMA10ExBin'] = roughECMA10ExBin
        dataByStim.loc[renderNames[ii], 'RoughECMA05ExBin'] = roughECMA05ExBin
        dataByStim.loc[renderNames[ii], 'RoughDW10ExMaxLR'] = roughDW10ExMaxLR
        dataByStim.loc[renderNames[ii], 'RoughDW05ExMaxLR'] = roughDW05ExMaxLR
        dataByStim.loc[renderNames[ii], 'FluctOV10ExMaxLR'] = fluctOV10ExMaxLR
        dataByStim.loc[renderNames[ii], 'FluctOV05ExMaxLR'] = fluctOV05ExMaxLR
        dataByStim.loc[renderNames[ii], 'SharpAurSHMPowAvgMaxLR'] = sharpASHMPowAvgMaxLR
        dataByStim.loc[renderNames[ii], 'SharpAurSHM05ExMaxLR'] = sharpASHM05ExMaxLR

        # calculation section for SQM differences
        # NOTE: THIS SECTION RELIES ON THE ALPHABETIC ORDER OF THE STIMULI FILES AS
        # ORIGINALLY NAMED: EACH AMBIENT SOUND FILE PRECEDING THE CORRESPONDING
        # COMBINED STIMULI FILES. IF THE ORDERING OR FILE NAMING IS CHANGED, THIS
        # CALCULATION WILL BE INVALIDATED AND *MAY ALSO* CAUSE AN ERROR.
        # a rolling 50 ms window is applied to average the SQM values over time -
        # this is to reduce uncertainty due to imperfect time-alignment between the
        # ambient vs combined stimuli files (all are from recordings, so there will
        # be some slippage due to imperfect editing)
        if renderNames[ii] in ["A1", "A2", "B2"]:
            # NOTE: we could dropna() the first <windowT values, but these will be
            # ignored anyway in the statistical analysis, assuming start_skipT >
            # windowT

            # calculate moving average values for ambient stimulus
            ambSpecTonalECMALMovAvg = specTonalECMAL.T.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean().T
            ambSpecTonalECMARMovAvg = specTonalECMAR.T.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean().T
            ambSpecTonLdECMALMovAvg = specTonLdECMAL.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()
            ambSpecTonLdECMARMovAvg = specTonLdECMAR.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()
            ambSpecRoughECMALMovAvg = specRoughECMAL.T.rolling(window=int(np.ceil(sampleRateRoughECMA*windowT))).mean().T
            ambSpecRoughECMARMovAvg = specRoughECMAR.T.rolling(window=int(np.ceil(sampleRateRoughECMA*windowT))).mean().T
            ambSpecFluctOVLMovAvg = specFluctOVL.T.rolling(window=int(np.ceil(sampleRateFluctOV*windowT))).mean().T
            ambSpecFluctOVRMovAvg = specFluctOVR.T.rolling(window=int(np.ceil(sampleRateFluctOV*windowT))).mean().T
            ambSharpASHMTDepMovAvg = sharpASHMTDep.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()

        elif renderNames[ii][0:3] in ["A1_", "A2_", "B2_"]:
            # calculate moving average values for combined stimulus
            specTonalECMALMovAvg = specTonalECMAL.T.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean().T
            specTonalECMARMovAvg = specTonalECMAR.T.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean().T
            specTonLdECMALMovAvg = specTonLdECMAL.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()
            specTonLdHMRSMovAvg = specTonLdECMAR.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()
            specRoughECMALMovAvg = specRoughECMAL.T.rolling(window=int(np.ceil(sampleRateRoughECMA*windowT))).mean().T
            specRoughECMARMovAvg = specRoughECMAR.T.rolling(window=int(np.ceil(sampleRateRoughECMA*windowT))).mean().T
            specFluctOVLMovAvg = specFluctOVL.T.rolling(window=int(np.ceil(sampleRateFluctOV*windowT))).mean().T
            specFluctOVRMovAvg = specFluctOVR.T.rolling(window=int(np.ceil(sampleRateFluctOV*windowT))).mean().T
            sharpASHMTDepMovAvg = sharpASHMTDep.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()

            # # calculate differences and make negative values 0
            dSpecTonalECMAL = np.maximum(specTonalECMALMovAvg
                                          - ambSpecTonalECMALMovAvg, 0)
            dSpecTonalECMAR = np.maximum(specTonalECMARMovAvg
                                          - ambSpecTonalECMARMovAvg, 0)
            dSpecTonLdECMAL = np.maximum(specTonLdECMALMovAvg
                                          - ambSpecTonLdECMALMovAvg, 0)
            dSpecTonLdECMAR = np.maximum(specTonLdHMRSMovAvg
                                          - ambSpecTonLdECMARMovAvg, 0)
            dSpecRoughECMAL = np.maximum(specRoughECMALMovAvg
                                          - ambSpecRoughECMALMovAvg, 0)
            dSpecRoughECMAR = np.maximum(specRoughECMARMovAvg
                                          - ambSpecRoughECMARMovAvg, 0)
            dSpecFluctOVL = np.maximum(specFluctOVLMovAvg
                                          - ambSpecFluctOVLMovAvg, 0)
            dSpecFluctOVR = np.maximum(specFluctOVRMovAvg
                                          - ambSpecFluctOVRMovAvg, 0)
            dSharpASHMTDep = np.maximum(sharpASHMTDepMovAvg
                                        - ambSharpASHMTDepMovAvg, 0)

            # calculate aggregated difference values

            # 2-channel time-dependent tonality (max, not integration)
            dTonalECMATDep = pd.concat([dSpecTonalECMAL.max(axis=0),
                                        dSpecTonalECMAR.max(axis=0)],
                                        axis=1)
            
            # 2-channel time-dependent integrated tonality
            dTonalSHMIntTDep = pd.concat([dSpecTonalECMAL.sum(axis=0)*bandDiff0p5*0.348088948583815,
                                          dSpecTonalECMAR.sum(axis=0)*bandDiff0p5*0.348088948583815],
                                          axis=1)

            # mask for start/end skip and values <= 0.02
            dTonalECMATDepMaskL = dTonalECMATDep.loc[(dTonalECMATDep.index.values
                                                      > start_skipT)
                                                      & (dTonalECMATDep.index.values
                                                        < dTonalECMATDep.index.values.max()
                                                        - end_skipT)
                                                      & (dTonalECMATDep.loc[:, 0].values
                                                        > 0.02), 0]
            dTonalECMATDepMaskR = dTonalECMATDep.loc[(dTonalECMATDep.index.values
                                                      > start_skipT)
                                                      & (dTonalECMATDep.index.values
                                                        < dTonalECMATDep.index.values.max()
                                                        - end_skipT)
                                                      & (dTonalECMATDep.loc[:, 1].values
                                                        > 0.02), 1]

            # mask for start/end skip and values <= 0.02
            # NOTE: uses tonality magnitude mask from ECMA tonality
            dTonalSHMIntTDepMaskL = dTonalSHMIntTDep.loc[(dTonalSHMIntTDep.index.values
                                                          > start_skipT)
                                                          & (dTonalSHMIntTDep.index.values
                                                            < dTonalSHMIntTDep.index.values.max()
                                                            - end_skipT)
                                                          & (dTonalECMATDep.loc[:, 0].values
                                                            > 0.02), 0]  # see NOTE above
            dTonalSHMIntTDepMaskR = dTonalSHMIntTDep.loc[(dTonalSHMIntTDep.index.values
                                                          > start_skipT)
                                                          & (dTonalSHMIntTDep.index.values
                                                            < dTonalSHMIntTDep.index.values.max()
                                                            - end_skipT)
                                                          & (dTonalECMATDep.loc[:, 1].values
                                                              > 0.02), 1]  # see NOTE above

            # 2-channel time-averaged tonality (omitting T<=0.02)
            dTonalECMAAvgL = dTonalECMATDepMaskL.mean(axis=0)
            dTonalECMAAvgR = dTonalECMATDepMaskL.mean(axis=0)
            dTonalECMAAvgMaxLR = max(dTonalECMAAvgL, dTonalECMAAvgR)
            # 2-channel 5% exceeded tonality (automatically omitting T<=0.02)
            dTonalECMA05ExL = dTonalECMATDepMaskL.quantile(q=0.95)
            dTonalECMA05ExR = dTonalECMATDepMaskL.quantile(q=0.95)
            dTonalECMA05ExMaxLR = max(dTonalECMA05ExL, dTonalECMA05ExR)

            # 2-channel time-averaged integrated tonality (omitting T<=0.02)
            dTonalSHMIntAvgL = dTonalSHMIntTDepMaskL.mean(axis=0)
            dTonalSHMIntAvgR = dTonalSHMIntTDepMaskL.mean(axis=0)
            dTonalSHMIntAvgMaxLR = max(dTonalSHMIntAvgL, dTonalSHMIntAvgR)
            # 2-channel 5% exceeded integated tonality (automatically omitting T<=0.02)
            dTonalSHMInt05ExL = dTonalSHMIntTDepMaskL.quantile(q=0.95)
            dTonalSHMInt05ExR = dTonalSHMIntTDepMaskL.quantile(q=0.95)
            dTonalSHMInt05ExMaxLR = max(dTonalSHMInt05ExL, dTonalSHMInt05ExR)

            # binaural specific tonal loudness (ECMA-418-2:2022 Equation 118)
            dSpecTonLdECMABin = ((dSpecTonLdECMAL**2
                                  + dSpecTonLdECMAR**2)/2).pow(0.5)
            # binaural time-dependent tonal loudness (ECMA-418-2:2022 Equation 116)
            dTonLdECMATDepBin = dSpecTonLdECMABin.sum(axis=0)*bandDiff0p5

            # mask for start/end skip
            dTonLdECMATDepBinMask = dTonLdECMATDepBin.loc[(dTonLdECMATDepBin.index.values
                                                            > start_skipT)
                                                          & (dTonLdECMATDepBin.index.values
                                                              < dTonLdECMATDepBin.index.values.max()
                                                              - end_skipT)]

            # binaural overall (power-averaged) tonal loudness (ECMA-418-2:2022
            # Equation 117)
            dTonLdECMAPowAvgBin = ((dTonLdECMATDepBinMask**(1/np.log10(2))).sum()
                                    / len(dTonLdECMATDepBinMask))**np.log10(2)

            # binaural 5% exceeded tonal loudness
            dTonLdECMA05ExBin = dTonLdECMATDepBinMask.quantile(q=0.95)

            # binaural specific roughness (ECMA-418-2:2022 Equation 112)
            dSpecRoughECMABin = ((dSpecRoughECMAL**2
                                  + dSpecRoughECMAL**2)/2).pow(0.5)

            # binaural time-dependent roughness
            dRoughECMATDepBin = dSpecRoughECMABin.sum(axis=0)*bandDiff0p5

            # mask for start/end skip
            dRoughECMATDepBinMask = dRoughECMATDepBin.loc[(dRoughECMATDepBin.index.values
                                                            > start_skipT)
                                                          & (dRoughECMATDepBin.index.values
                                                              < dRoughECMATDepBin.index.values.max()
                                                              - end_skipT)]

            # binaural overall (90th percentile = 10% exceeded) roughness
            dRoughECMA10ExBin = dRoughECMATDepBinMask.quantile(q=0.90)
            # binaural overall (95th percentile = 5% exceeded) roughness
            dRoughECMA05ExBin = dRoughECMATDepBinMask.quantile(q=0.95)
            
            # 2-channel time-dependent integrated fluctuation strength
            dFluctOVTDep = pd.concat([dSpecFluctOVL.sum(axis=0)*bandDiff0p5,
                                      dSpecFluctOVR.sum(axis=0)*bandDiff0p5],
                                      axis=1)
            
            # 2-channel fluctuation strength masked for start/end skip
            dFluctOVTDepMask = dFluctOVTDep.loc[(dFluctOVTDep.index.values
                                                  > start_skipT).transpose()
                                                & (dFluctOVTDep.index.values
                                                    < dFluctOVTDep.index.values.max()
                                                    - end_skipT).transpose()]
            
            # 2-channel overall 5% exceeded fluctuation strength 
            dFluctOV05Ex = dFluctOVTDepMask.quantile(q=0.95)
            # max of l/r channel overall 5% exceeded fluctuation strength 
            dFluctOV05ExMaxLR = dFluctOV05Ex.max()
            # 2-channel overall 10% exceeded fluctuation strength 
            dFluctOV10Ex = dFluctOVTDepMask.quantile(q=0.90)
            # max of l/r channel overall 10% exceeded fluctuation strength 
            dFluctOV10ExMaxLR = dFluctOV10Ex.max()

            # 2-channel sharpness masked for start/end skip
            dSharpASHMTDepMask = dSharpASHMTDep.loc[(dSharpASHMTDep.index.values
                                                      > start_skipT).transpose()
                                                    & (dSharpASHMTDep.index.values
                                                        < dSharpASHMTDep.index.values.max()
                                                        - end_skipT).transpose()]

            # 2-channel overall (power-averaged) sharpness
            dSharpASHMPowAvg = ((dSharpASHMTDepMask**(1/np.log10(2))).sum(axis=0)
                                / len(dSharpASHMTDepMask))**np.log10(2)
            # max of l/r channel overall (power-averaged) sharpness
            dSharpASHMPowAvgMaxLR = dSharpASHMPowAvg.max()
            # 2-channel overall 5% exceeded sharpness
            dSharpASHM05Ex = dSharpASHMTDepMask.quantile(q=0.95)
            # max of l/r channel overall 5% exceeded sharpness
            dSharpASHM05ExMaxLR = dSharpASHM05Ex.max()

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'dTonalECMAAvgMaxLR'] = dTonalECMAAvgMaxLR
            dataByStim.loc[renderNames[ii], 'dTonalECMA05ExMaxLR'] = dTonalECMA05ExMaxLR
            dataByStim.loc[renderNames[ii], 'dTonalSHMIntAvgMaxLR'] = dTonalSHMIntAvgMaxLR
            dataByStim.loc[renderNames[ii], 'dTonalSHMInt05ExMaxLR'] = dTonalSHMInt05ExMaxLR
            dataByStim.loc[renderNames[ii], 'dTonLdECMAPowAvgBin'] = dTonLdECMAPowAvgBin
            dataByStim.loc[renderNames[ii], 'dTonLdECMA05ExBin'] = dTonLdECMA05ExBin
            dataByStim.loc[renderNames[ii], 'dRoughECMA10ExBin'] = dRoughECMA10ExBin
            dataByStim.loc[renderNames[ii], 'dRoughECMA05ExBin'] = dRoughECMA05ExBin
            dataByStim.loc[renderNames[ii], 'dFluctOV10ExMaxLR'] = dFluctOV10ExMaxLR
            dataByStim.loc[renderNames[ii], 'dFluctOV05ExMaxLR'] = dFluctOV05ExMaxLR
            dataByStim.loc[renderNames[ii], 'dSharpAurSHMPowAvgMaxLR'] = dSharpASHMPowAvgMaxLR
            dataByStim.loc[renderNames[ii], 'dSharpAurSHM05ExMaxLR'] = dSharpASHM05ExMaxLR

# end of for loop over MATLAB SQM files

# From ArtemiS calculations
# -------------------------

# open xlsx file selection dialog and assign filepaths to list
# PROJECT NOTE: the calculated files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Analysis/ArtemiS/CALHEQ?csf=1&web=1&e=R2uFcs
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance()

fileExtsXLSX = "*.xlsx"
filelistXLSX = list(QFileDialog.getOpenFileNames(filter=fileExtsXLSX,
                                                 caption=r"Select ArtemiS output files in '03 Experiment\Experiment 1\Analysis\ArtemiS\CALHEQ\xlsx'"))[0]
filelistXLSX.sort()

fileExtsASC = "*.asc"
filelistASC = list(QFileDialog.getOpenFileNames(filter=fileExtsASC,
                                                caption=r"Select ArtemiS output files in '03 Experiment\Experiment 1\Analysis\ArtemiS\CALHEQ\asc'"))[0]

filelistASC.sort()

filelist = filelistXLSX + filelistASC

filenames = [filepath.split('/')[-1] for filepath in filelist]
renderNames = [filename.split('_CALHEQ_Pa')[0] for filename in filenames]

# loop over files to analyse
for ii, file in enumerate(filelist):
    if ii == 0:
        print("Processing sound quality metrics...\n")
    print(file.split('/')[-1])

    if file[-3:] == "asc":

        if file.find("Specific") != -1:
            filedata = np.array(pd.read_csv(file, sep="	", header=None))
            chan2row = np.arange(0, filedata.shape[0])[np.isnan(filedata[:, 0])][-1]
    
            # Calculate overall Sottek Hearing Model impulsiveness from specific
            # time-dependent impulsiveness
            specImpulsSHML = pd.DataFrame(filedata[1:chan2row, 1:],
                                          index=filedata[1:chan2row, 0],
                                          columns=filedata[0, 1:])
            specImpulsSHMR = pd.DataFrame(filedata[chan2row + 1:, 1:],
                                          index=filedata[chan2row + 1:, 0],
                                          columns=filedata[chan2row, 1:])
    
            # calculation section for SQM differences
            # NOTE: THIS SECTION RELIES ON THE ALPHABETIC ORDER OF THE STIMULI FILES AS
            # ORIGINALLY NAMED: EACH AMBIENT SOUND FILE PRECEDING THE CORRESPONDING
            # COMBINED STIMULI FILES. IF THE ORDERING OR FILE NAMING IS CHANGED, THIS
            # CALCULATION WILL BE INVALIDATED AND *MAY ALSO* CAUSE AN ERROR.
            # a rolling 50 ms window is applied to average the SQM values over time -
            # this is to reduce uncertainty due to imperfect time-alignment between the
            # ambient vs combined stimuli files (all are from recordings, so there will
            # be some slippage due to imperfect editing)
            if renderNames[ii] in ["A1", "A2", "B2"]:
                # NOTE: we could dropna() the first <windowT values, but these will be
                # ignored anyway in the statistical analysis, assuming start_skipT >
                # windowT
    
                # calculate moving average values for ambient stimulus
                ambSpecImpulsSHMLMovAvg = specImpulsSHML.rolling(window=int(np.ceil(sampleRateImpulsSHM*windowT))).mean().T
                ambSpecImpulsSHMRMovAvg = specImpulsSHMR.rolling(window=int(np.ceil(sampleRateImpulsSHM*windowT))).mean().T
    
            elif renderNames[ii][0:3] in ["A1_", "A2_", "B2_"]:
                # calculate moving average values for combined stimulus
                specImpulsSHMLMovAvg = specImpulsSHML.rolling(window=int(np.ceil(sampleRateImpulsSHM*windowT))).mean().T
                specImpulsSHMRMovAvg = specImpulsSHMR.rolling(window=int(np.ceil(sampleRateImpulsSHM*windowT))).mean().T
    
                # calculate differences and make negative values 0
                dImpulsSHMTDep = pd.concat([np.maximum(specImpulsSHMLMovAvg
                                           - ambSpecImpulsSHMLMovAvg, 0).sum(axis=0).to_frame(name="Channel 1"),
                                            np.maximum(specImpulsSHMRMovAvg
                                           - ambSpecImpulsSHMRMovAvg, 0).sum(axis=0).to_frame(name="Channel 2")],
                                           axis=1)
    
                # calculate aggregated difference values
                # 2-channel impulsiveness masked for start/end skip
                dImpulsSHMTDepMask = dImpulsSHMTDep.loc[(dImpulsSHMTDep.index.values
                                                         > start_skipT).transpose()
                                                        & (dImpulsSHMTDep.index.values
                                                           < dImpulsSHMTDep.index.values.max()
                                                           - end_skipT).transpose()]
                # 2-channel overall (power-averaged) impulsiveness
                dImpulsSHMPowAvg = ((dImpulsSHMTDepMask**(1/np.log10(2))).sum(axis=0)
                                    / len(dImpulsSHMTDepMask))**np.log10(2)
                # max of l/r overall (power-averaged) impulsiveness
                dImpulsSHMPowAvgMaxLR = dImpulsSHMPowAvg.max()
                # 2-channel overall (mean) impulsiveness
                dImpulsSHMMean = dImpulsSHMTDepMask.mean(axis=0)
                # max of l/r channel overall (mean) impulsiveness
                dImpulsSHMMeanMaxLR = dImpulsSHMMean.max()
    
                # 2-channel overall 5% exceeded impulsiveness
                dImpulsSHM05Ex = dImpulsSHMTDepMask.quantile(q=0.95)
                # max of l/r channel overall 5% exceeded impulsiveness
                dImpulsSHM05ExMaxLR = dImpulsSHM05Ex.max()
    
                # add results to output DataFrame
                dataByStim.loc[renderNames[ii], 'dImpulsSHMPowAvgMaxLR'] = dImpulsSHMPowAvgMaxLR
                dataByStim.loc[renderNames[ii], 'dImpulsSHMAvgMaxLR'] = dImpulsSHMMeanMaxLR
                dataByStim.loc[renderNames[ii], 'dImpulsSHM05ExMaxLR'] = dImpulsSHM05ExMaxLR
            
        # if-else branch for not specific impulsiveness (overall)
        else:
            
            # read in time-dependent impulsiveness
            impulsSHMTDep = pd.read_csv(file, sep="	", header=None, index_col=0)
            
            # mask for start/end skip and 0 values
            impulsSHMTDepMask = impulsSHMTDep.loc[(impulsSHMTDep.index.values
                                                   > start_skipT).transpose()
                                                  & (impulsSHMTDep.index.values
                                                     < impulsSHMTDep.index.values.max()
                                                     - end_skipT).transpose()]

            # 2-channel overall (power-averaged) impulsiveness
            impulsSHMPowAvg = ((impulsSHMTDepMask**(1/np.log10(2))).sum(axis=0)
                               / len(impulsSHMTDepMask))**np.log10(2)
            # max of l/r overall (power-averaged) impulsiveness
            impulsSHMPowAvgMaxLR = impulsSHMPowAvg.max()
            # 2-channel overall (mean) impulsiveness
            impulsSHMMean = impulsSHMTDepMask.mean(axis=0)
            # max of l/r channel overall (mean) impulsiveness
            impulsSHMMeanMaxLR = impulsSHMMean.max()

            # 2-channel overall 5% exceeded impulsiveness
            impulsSHM05Ex = impulsSHMTDepMask.quantile(q=0.95)

            # max of l/r channel overall 5% exceeded impulsiveness
            impulsSHM05ExMaxLR = impulsSHM05Ex.max()

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'ImpulsSHMPowAvgMaxLR'] = impulsSHMPowAvgMaxLR
            dataByStim.loc[renderNames[ii], 'ImpulsSHMAvgMaxLR'] = impulsSHMMeanMaxLR
            dataByStim.loc[renderNames[ii], 'ImpulsSHM05ExMaxLR'] = impulsSHM05ExMaxLR



    # end of if branch for .asc files

    elif file[-4:] == "xlsx":
        workbookdata = pd.read_excel(io=file, sheet_name=None)

        recFile = workbookdata['Sheet1'].columns[5].split(sep='\'')[1]

        # check file is correct and added sheets match first sheet
        if filenames[ii].split(".xlsx")[0] != recFile.split(".wav")[0]:
            raise ValueError("Error: The recording analysis does not match the filename. Check file " + filenames[ii])
        elif recFile != workbookdata['Sheet8'].columns[5].split(sep='\'')[1]:
            raise ValueError("Error: The recording analysis does not match across sheets. Check file " + filenames[ii])
        elif recFile != workbookdata['Sheet9'].columns[5].split(sep='\'')[1]:
            raise ValueError("Error: The recording analysis does not match across sheets. Check file " + filenames[ii])

        # fix any problem filenames (legacy: should be unnecessary)
        if recFile.find("YnTy") != -1:
            recFile = recFile.replace("YnTy", "H520")
        if recFile.find(".wav_Pa.wav") != -1:
            recFile = recFile.replace(".wav_Pa.wav", "_Pa.wav")

        # Calculate overall Sottek Hearing Model fluctuation strength from
        # 2-channel specific fluctuation strength
        # calculate according to ECMA-418-2:2022 approach for roughness
        # left channel
        specFluctSHML = pd.DataFrame(workbookdata['Sheet6'].iloc[14:, 1:].values,
                                     columns=workbookdata['Sheet6'].iloc[13, 1:],
                                     index=workbookdata['Sheet6'].iloc[14:, 0])
        # right channel
        specFluctSHMR = pd.DataFrame(workbookdata['Sheet7'].iloc[14:, 1:].values,
                                     columns=workbookdata['Sheet7'].iloc[13, 1:],
                                     index=workbookdata['Sheet7'].iloc[14:, 0])

        # binaural specific fluctuation strength
        # (using ECMA-418-2:2022 Equation 112 for roughness)
        specFluctSHMBin = ((specFluctSHML**2
                            + specFluctSHMR**2)/2).pow(0.5)
        # binaural time-dependent fluctuation strength
        # (using ECMA-418-2:2022 Equation 111 for roughness)
        fluctSHMTDepBin = specFluctSHMBin.sum(axis=0)*bandDiff0p5

        # mask for start/end skip
        fluctSHMTDepBinMask = fluctSHMTDepBin.loc[(fluctSHMTDepBin.index.values
                                                   > start_skipT)
                                                  & (fluctSHMTDepBin.index.values
                                                     < fluctSHMTDepBin.index.values.max()
                                                     - end_skipT)]

        # binaural overall (90th percentile = 10% exceeded) fluctuation strength
        # (using ECMA-418-2:2022 Section 7.1.8 for roughness)
        fluctSHMTime10ExBin = fluctSHMTDepBinMask.quantile(q=0.90)
        # binaural overall (95th percentile = 5% exceeded) fluctuation strength
        fluctSHMTime05ExBin = fluctSHMTDepBinMask.quantile(q=0.95)

        # Calculate overall ISO 532-1 loudness from 2-channel time-varing loudness
        loudISO1TDep = pd.DataFrame(workbookdata['Sheet1'].iloc[13:, 1:3].values,
                                    columns=workbookdata['Sheet1'].iloc[12, 1:3],
                                    index=workbookdata['Sheet1'].iloc[13:, 0])

        # mask for start/end skip and 0 values
        loudISO1TDepMask = loudISO1TDep.loc[(loudISO1TDep.index.values
                                             > start_skipT).transpose()
                                            & (loudISO1TDep.index.values
                                               < loudISO1TDep.index.values.max()
                                               - end_skipT).transpose()]

        # 2-channel overall (5% exceeded = 95th percentile) loudness
        loudISO105Ex = loudISO1TDepMask.quantile(q=0.95)
        # max of l/r channel (5% exceeded = 95th percentile) loudness
        loudISO105ExMaxLR = loudISO105Ex.max()

        # 2-channel overall (power-averaged) loudness
        loudISO1PowAvg = ((loudISO1TDepMask**(1/np.log10(2))).sum(axis=0)
                          / len(loudISO1TDepMask))**np.log10(2)
        # max of l/r channel (95th-percentile) loudness
        loudISO1PowAvgMaxLR = loudISO1PowAvg.max()

        # Calculate overall ISO 532-3 loudness from binaural time-varing loudness
        loudISO3TDepBin = pd.DataFrame(workbookdata['Sheet9'].iloc[13:, 1:2].values,
                                       columns=workbookdata['Sheet9'].iloc[12, 1:2],
                                       index=workbookdata['Sheet9'].iloc[13:, 0])

        # mask for start/end skip
        loudISO3TDepBinMask = loudISO3TDepBin.loc[(loudISO3TDepBin.index.values
                                                   > start_skipT)
                                                  & (loudISO3TDepBin.index.values
                                                     < loudISO3TDepBin.index.values.max()
                                                     - end_skipT)]

        # binaural overall (power-averaged) loudness
        loudISO3PowAvgBin = ((loudISO3TDepBinMask**(1/np.log10(2))).sum(axis=0)
                             / len(loudISO3TDepBinMask))**np.log10(2)
        loudISO3PowAvgBin = loudISO3PowAvgBin.iloc[0]  # convert series to float

        # Calculate overall DIN 45692 sharpness from 2-channel time-dependent
        # sharpness
        sharpDINTDep = pd.DataFrame(workbookdata['Sheet3'].iloc[13:, 1:3].values,
                                    columns=workbookdata['Sheet3'].iloc[12, 1:3],
                                    index=workbookdata['Sheet3'].iloc[13:, 0])

        # mask for start/end skip
        sharpDINTDepMask = sharpDINTDep.loc[(sharpDINTDep.index.values
                                             > start_skipT)
                                            & (sharpDINTDep.index.values
                                               < sharpDINTDep.index.values.max()
                                               - end_skipT), :]

        # 2-channel overall (power-averaged) sharpness
        sharpDINPowAvg = ((sharpDINTDepMask**(1/np.log10(2))).sum(axis=0)
                          / len(sharpDINTDepMask))**np.log10(2)
        # max of l/r channel overall (power-averaged) sharpness
        sharpDINPowAvgMaxLR = sharpDINPowAvg.max()
        # 2-channel overall 5% exceeded sharpness
        sharpDIN05Ex = sharpDINTDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded sharpness
        sharpDIN05ExMaxLR = sharpDIN05Ex.max()

        # Calculate overall von Bismarck | ISO 532-1 sharpness from 2-channel
        # time-dependent sharpness
        sharpvBISO1TDep = pd.DataFrame(workbookdata['Sheet4'].iloc[13:, 1:3].values,
                                       columns=workbookdata['Sheet4'].iloc[12, 1:3],
                                       index=workbookdata['Sheet4'].iloc[13:, 0])

        # mask for start/end skip
        sharpvBISO1TDepMask = sharpvBISO1TDep.loc[(sharpvBISO1TDep.index.values
                                                   > start_skipT)
                                                  & (sharpvBISO1TDep.index.values
                                                     < sharpvBISO1TDep.index.values.max()
                                                     - end_skipT), :]

        # 2-channel overall (power-averaged) sharpness
        sharpvBISO1PowAvg = ((sharpvBISO1TDepMask**(1/np.log10(2))).sum(axis=0)
                             / len(sharpvBISO1TDepMask))**np.log10(2)
        # max of l/r channel overall (power-averaged) sharpness
        sharpvBISO1PowAvgMaxLR = sharpvBISO1PowAvg.max()
        # 2-channel overall 5% exceeded sharpness
        sharpvBISO105Ex = sharpvBISO1TDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded sharpness
        sharpvBISO105ExMaxLR = sharpvBISO105Ex.max()

        # Calculate overall Aures+ISO532-3 sharpness from 2-channel time-dependent
        # sharpness
        sharpAISO3TDep = pd.DataFrame(workbookdata['Sheet8'].iloc[13:, 1:3].values,
                                      columns=workbookdata['Sheet8'].iloc[12, 1:3],
                                      index=workbookdata['Sheet8'].iloc[13:, 0])

        # mask for start/end skip
        sharpAISO3TDepMask = sharpAISO3TDep.loc[(sharpAISO3TDep.index.values
                                                 > start_skipT)
                                                & (sharpAISO3TDep.index.values
                                                   < sharpAISO3TDep.index.values.max()
                                                   - end_skipT), :]

        # 2-channel overall (power-averaged) sharpness
        sharpAISO3PowAvg = ((sharpAISO3TDepMask**(1/np.log10(2))).sum(axis=0)
                            / len(sharpAISO3TDepMask))**np.log10(2)
        # max of l/r channel overall (power-averaged) sharpness
        sharpAISO3PowAvgMaxLR = sharpAISO3PowAvg.max()
        # 2-channel overall 5% exceeded sharpness
        sharpAISO305Ex = sharpAISO3TDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded sharpness
        sharpAISO305ExMaxLR = sharpAISO305Ex.max()

        # Calculate overall Aures+ISO532-1 sharpness from 2-channel time-dependent
        # sharpness
        sharpAISO1TDep = pd.DataFrame(workbookdata['Sheet5'].iloc[13:, 1:3].values,
                                      columns=workbookdata['Sheet5'].iloc[12, 1:3],
                                      index=workbookdata['Sheet5'].iloc[13:, 0])

        # mask for start/end skip
        sharpAISO1TDepMask = sharpAISO1TDep.loc[(sharpAISO1TDep.index.values
                                                 > start_skipT)
                                                & (sharpAISO1TDep.index.values
                                                   < sharpAISO1TDep.index.values.max()
                                                   - end_skipT), :]

        # 2-channel overall (power-averaged) sharpness
        sharpAISO1PowAvg = ((sharpAISO1TDepMask**(1/np.log10(2))).sum(axis=0)
                            / len(sharpAISO1TDepMask))**np.log10(2)
        # max of l/r channel overall (power-averaged) sharpness
        sharpAISO1PowAvgMaxLR = sharpAISO1PowAvg.max()
        # 2-channel overall 5% exceeded sharpness
        sharpAISO105Ex = sharpAISO1TDepMask.quantile(q=0.95)
        # max of l/r channel overall 5% exceeded sharpness
        sharpAISO105ExMaxLR = sharpAISO105Ex.max()

        # add results to output DataFrame
        dataByStim.loc[renderNames[ii], 'FluctSHM10ExBin'] = fluctSHMTime10ExBin
        dataByStim.loc[renderNames[ii], 'FluctSHM05ExBin'] = fluctSHMTime05ExBin
        dataByStim.loc[renderNames[ii], 'LoudISO105ExMaxLR'] = loudISO105ExMaxLR
        dataByStim.loc[renderNames[ii], 'LoudISO1PowAvgMaxLR'] = loudISO1PowAvgMaxLR
        dataByStim.loc[renderNames[ii], 'LoudISO3PowAvgBin'] = loudISO3PowAvgBin
        dataByStim.loc[renderNames[ii], 'SharpAurISO3PowAvgMaxLR'] = sharpAISO3PowAvgMaxLR
        dataByStim.loc[renderNames[ii], 'SharpAurISO305ExMaxLR'] = sharpAISO305ExMaxLR
        dataByStim.loc[renderNames[ii], 'SharpAurISO1PowAvgMaxLR'] = sharpAISO1PowAvgMaxLR
        dataByStim.loc[renderNames[ii], 'SharpAurISO105ExMaxLR'] = sharpAISO105ExMaxLR
        dataByStim.loc[renderNames[ii], 'SharpvBISO1PowAvgMaxLR'] = sharpvBISO1PowAvgMaxLR
        dataByStim.loc[renderNames[ii], 'SharpvBISO105ExMaxLR'] = sharpvBISO105ExMaxLR
        dataByStim.loc[renderNames[ii], 'SharpDINPowAvgMaxLR'] = sharpDINPowAvgMaxLR
        dataByStim.loc[renderNames[ii], 'SharpDIN05ExMaxLR'] = sharpDIN05ExMaxLR

        # calculation section for SQM differences
        # NOTE: THIS SECTION RELIES ON THE ALPHABETIC ORDER OF THE STIMULI FILES AS
        # ORIGINALLY NAMED: EACH AMBIENT SOUND FILE PRECEDING THE CORRESPONDING
        # COMBINED STIMULI FILES. IF THE ORDERING OR FILE NAMING IS CHANGED, THIS
        # CALCULATION WILL BE INVALIDATED AND *MAY ALSO* CAUSE AN ERROR.
        # a rolling 50 ms window is applied to average the SQM values over time -
        # this is to reduce uncertainty due to imperfect time-alignment between the
        # ambient vs combined stimuli files (all are from recordings, so there will
        # be some slippage due to imperfect editing)
        if renderNames[ii] in ["A1", "A2", "B2"]:
            # NOTE: we could dropna() the first <windowT values, but these will be
            # ignored anyway in the statistical analysis, assuming start_skipT >
            # windowT

            # calculate moving average values for ambient stimulus
            ambSpecFluctSHMLMovAvg = specFluctSHML.T.rolling(window=int(np.ceil((sampleRateFluctSHM)*windowT))).mean().T
            ambSpecFluctSHMRMovAvg = specFluctSHMR.T.rolling(window=int(np.ceil((sampleRateFluctSHM)*windowT))).mean().T
            ambSharpAISO3TDepMovAvg = sharpAISO3TDep.rolling(window=int(np.ceil(sampleRateLoudISO3*windowT))).mean()

        elif renderNames[ii][0:3] in ["A1_", "A2_", "B2_"]:
            # calculate moving average values for combined stimulus
            specFluctSHMLMovAvg = specFluctSHML.T.rolling(window=int(np.ceil((sampleRateFluctSHM)*windowT))).mean().T
            specFluctSHMRMovAvg = specFluctSHMR.T.rolling(window=int(np.ceil((sampleRateFluctSHM)*windowT))).mean().T
            sharpAISO3TDepMovAvg = sharpAISO3TDep.rolling(window=int(np.ceil(sampleRateLoudISO3*windowT))).mean()

            # calculate differences and make negative values 0
            dSpecFluctSHML = np.maximum(specFluctSHMLMovAvg
                                        - ambSpecFluctSHMLMovAvg, 0)
            dSpecFluctSHMR = np.maximum(specFluctSHMRMovAvg
                                        - ambSpecFluctSHMRMovAvg, 0)
            dSharpAISO3TDep = np.maximum(sharpAISO3TDepMovAvg
                                         - ambSharpAISO3TDepMovAvg, 0)

            # calculate aggregated difference values

            # binaural specific fluctuation strength
            # (using ECMA-418-2:2022 Equation 112 for roughness)
            dSpecFluctSHMBin = ((dSpecFluctSHML**2
                                 + dSpecFluctSHMR**2)/2).pow(0.2)
            # binaural time-dependent fluctuation strength
            # (using ECMA-418-2:2022 Equation 111 for roughness)
            dFluctSHMTDepBin = dSpecFluctSHMBin.sum(axis=0)*bandDiff0p5

            # mask for start/end skip
            dFluctSHMTDepBinMask = dFluctSHMTDepBin.loc[(dFluctSHMTDepBin.index.values
                                                         > start_skipT)
                                                        & (dFluctSHMTDepBin.index.values
                                                           < dFluctSHMTDepBin.index.values.max()
                                                           - end_skipT)]

            # binaural overall (90th percentile = 10% exceeded) fluctuation strength
            # (using ECMA-418-2:2022 Section 7.1.8 for roughness)
            dFluctSHM10ExBin = dFluctSHMTDepBinMask.quantile(q=0.90)
            # binaural overall (95th percentile = 5% exceeded) fluctuation strength
            dFluctSHM05ExBin = dFluctSHMTDepBinMask.quantile(q=0.95)

            # 2-channel sharpness masked for start/end skip
            dSharpAISO3TDepMask = dSharpAISO3TDep.loc[(dSharpAISO3TDep.index.values
                                                       > start_skipT).transpose()
                                                      & (dSharpAISO3TDep.index.values
                                                         < dSharpAISO3TDep.index.values.max()
                                                         - end_skipT).transpose()]

            # 2-channel overall (power-averaged) sharpness
            dSharpAISO3PowAvg = ((dSharpAISO3TDepMask**(1/np.log10(2))).sum(axis=0)
                                 / len(dSharpAISO3TDepMask))**np.log10(2)
            # max of l/r channel overall (power-averaged) sharpness
            dSharpAISO3PowAvgMaxLR = dSharpAISO3PowAvg.max()
            # 2-channel overall 5% exceeded sharpness
            dSharpAISO305Ex = dSharpAISO3TDepMask.quantile(q=0.95)
            # max of l/r channel overall 5% exceeded sharpness
            dSharpAISO305ExMaxLR = dSharpAISO305Ex.max()

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'dFluctSHM10ExBin'] = dFluctSHM10ExBin
            dataByStim.loc[renderNames[ii], 'dFluctSHM05ExBin'] = dFluctSHM05ExBin
            dataByStim.loc[renderNames[ii], 'dSharpAurISO3PowAvgMaxLR'] = dSharpAISO3PowAvgMaxLR
            dataByStim.loc[renderNames[ii], 'dSharpAurISO305ExMaxLR'] = dSharpAISO305ExMaxLR

# end of for loop over ArtemiS SQM files

# rearrange SQM columns
dataByStim[indicesPsycho] = dataByStim.reindex(columns=indicesPsycho)

# add UAS-only and ambient-only data to combined stimuli
# ------------------------------------------------------
indicesDiffPsycho = ["dTonalECMAAvgMaxLR",
                     "dTonalECMA05ExMaxLR",
                     "dTonalSHMIntAvgMaxLR",
                     "dTonalSHMInt05ExMaxLR",
                     "dTonLdECMAPowAvgBin",
                     "dTonLdECMA05ExBin",
                     "dRoughECMA10ExBin",
                     "dRoughECMA05ExBin",
                     "dRoughFZ10ExMaxLR",
                     "dRoughFZ05ExMaxLR",
                     "dFluctSHM10ExBin",
                     "dFluctSHM05ExBin",
                     "dFluctOV10ExMaxLR",
                     "dFluctOV05ExMaxLR",
                     "dSharpAurSHMPowAvgMaxLR",
                     "dSharpAurSHM05ExMaxLR",
                     "dSharpAurISO3PowAvgMaxLR",
                     "dSharpAurISO305ExMaxLR",
                     "dImpulsSHMPowAvgMaxLR",
                     "dImpulsSHMAvgMaxLR",
                     "dImpulsSHM05ExMaxLR"]

indicesAbsPsycho = [index for index in indicesPsycho
                    if index not in indicesDiffPsycho]

indices = indicesAcoustic + indicesPNL + indicesAbsPsycho

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
mask = dataByStim.CALBINRecFiles.str.find("A1_CALBIN") != -1
AmbonlyPtA1 = dataByStim.loc[mask, indices].copy()
AmbonlyPtA1 = AmbonlyPtA1.add_prefix("Amb", axis=1)
AmbonlyPtA1 = AmbonlyPtA1.loc[AmbonlyPtA1.index.repeat(
                              sum(dataByStim.index.str.find("A_") != -1) + 1)]
AmbonlyPtA1['newindex'] = UASonlyPtA1.index.union(dataByStim.index[dataByStim.CALBINRecFiles.str.find("A1_CALBIN") != -1])
AmbonlyPtA1.set_index('newindex', inplace=True)

# Part A2 ambient only
mask = dataByStim.CALBINRecFiles.str.find("A2_CALBIN") != -1
AmbonlyPtA2 = dataByStim.loc[mask, indices].copy()
AmbonlyPtA2 = AmbonlyPtA2.add_prefix("Amb", axis=1)
AmbonlyPtA2 = AmbonlyPtA2.loc[AmbonlyPtA2.index.repeat(
                              sum(dataByStim.index.str.find("A_") != -1) + 1)]
AmbonlyPtA2['newindex'] = UASonlyPtA2.index.union(dataByStim.index[dataByStim.CALBINRecFiles.str.find("A2_CALBIN") != -1])
AmbonlyPtA2.set_index('newindex', inplace=True)

AmbonlyPtA = pd.concat([AmbonlyPtA1, AmbonlyPtA2], axis=0)

# Part B ambient only
mask = dataByStim.CALBINRecFiles.str.find("B2_CALBIN") != -1
AmbonlyPtB = dataByStim.loc[mask, indices].copy()
AmbonlyPtB = AmbonlyPtB.add_prefix("Amb", axis=1)
AmbonlyPtB = AmbonlyPtB.loc[AmbonlyPtB.index.repeat(
                            sum(dataByStim.index.str.find("B_") != -1) + 1)]
AmbonlyPtB['newindex'] = UASonlyPtB.index.union(dataByStim.index[dataByStim.CALBINRecFiles.str.find("B2_CALBIN") != -1])
AmbonlyPtB.set_index('newindex', inplace=True)

# concatenate UAS only
Ambonly = pd.concat([AmbonlyPtA, AmbonlyPtB], axis=0)

# concatenate UAS only and ambient only
UASAmb = pd.concat([UASonly, Ambonly], axis=1)

# insert zeros for No UAS stimuli UAS SQMs
NoUASStims = ['A1', 'A2', 'B2']
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
dataByStim['LAELAF90diff'] = dataByStim['UASLAEMaxLR'] - dataByStim['AmbLAF90ExMaxLR']
dataByStim['LAELAF50diff'] = dataByStim['UASLAEMaxLR'] - dataByStim['AmbLAF50ExMaxLR']
dataByStim['LAELAeqdiff'] = dataByStim['UASLAEMaxLR'] - dataByStim['AmbLAeqMaxLR']
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
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%201/Analysis/deeuu_loudness/Partial_loudness?csf=1&web=1&e=fVScn8
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

fileExts = "*TermPartialLoudness.pkl"
filelist = list(QFileDialog.getOpenFileNames(filter=fileExts,
                                             caption="Open partial loudness files in '03 Experiment\Experiment 1\Analysis\deeuu_loudness\Partial_loudness'"))[0]
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

# reindex partial loudness DataFrame to match dataByStim index and merge into
# output DataFrame
partLoudnessST['newindex'] = [file.replace("_CALHEQ_PaShortTermPartialLoudness.pkl", "")
                              for file in list(partLoudnessST.index)]
partLoudnessST.set_index('newindex', inplace=True)
partLoudnessLT['newindex'] = [file.replace("_CALHEQ_PaLongTermPartialLoudness.pkl", "")
                              for file in list(partLoudnessLT.index)]
partLoudnessLT.set_index('newindex', inplace=True)

partLoudness = partLoudnessST.merge(right=partLoudnessLT, how='outer',
                                    left_index=True, right_index=True)

dataByStim = dataByStim.merge(partLoudness.astype(float), how='outer',
                              left_index=True, right_index=True)

# insert zeros for No UAS stimuli partial loudness
dataByStim.loc[NoUASStims, ['PartLoudMGSTPowAvg', 'PartLoudMGST05Ex',
                            'PartLoudMGLTPowAvg', 'PartLoudMGLT05Ex'] ] = 0

# --------------------------------------------------
# Detection-discounted UAS sound levels LAeq and LAE
# --------------------------------------------------

dataByStim['UASDisc0p5LAeqMaxLR'] = dataByStim['UASLAeqMaxLR'] - dataByStim['Detect0p5dBADiscMaxLR']
dataByStim['UASDisc0p5LAEMaxLR'] = dataByStim['UASLAEMaxLR'] - dataByStim['Detect0p5dBADiscMaxLR']
dataByStim['UASDisc0p1LAeqMaxLR'] = dataByStim['UASLAeqMaxLR'] - dataByStim['Detect0p1dBADiscMaxLR']
dataByStim['UASDisc0p1LAEMaxLR'] = dataByStim['UASLAEMaxLR'] - dataByStim['Detect0p1dBADiscMaxLR']

indicesDetect = indicesDetect + ['UASDisc0p5LAeqMaxLR', 'UASDisc0p5LAEMaxLR',
                                 'UASDisc0p1LAeqMaxLR', 'UASDisc0p1LAEMaxLR']


# reorganise column order

# move detectability metrics to after partial loudness metrics
indicesDetect.reverse()
for col in indicesDetect:
    popcol = dataByStim.pop(col)
    dataByStim.insert(loc=(dataByStim.columns.get_loc(partLoudness.columns[-1])
                           + 1),
                      column=col, value=popcol)
indicesDetect.reverse()  # revert list for use later   


# move SQM difference metrics to after detectability metrics
indicesDiffPsycho.reverse()
for col in indicesDiffPsycho:
    popcol = dataByStim.pop(col)
    dataByStim.insert(loc=(dataByStim.columns.get_loc(indicesDetect[-1])
                           + 1),
                      column=col, value=popcol)
indicesDiffPsycho.reverse()  # revert list for use later   



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
                                       caption=r"Select test response data files in '03 Experiment\Experiment 1\Test_files\Results'")[0]

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
        partAResponses.loc[(partAResponses['ID#'] == iD)
                           & (partAResponses['Stim File'].str.contains("A1")),
                           "d" + response] = (partAResponses.loc[(partAResponses['ID#']
                                                                  == iD)
                                                                 & (partAResponses['Stim File'].str.contains("A1"))].loc[:,
                                                                                                                         response]
                                              - partAResponses.loc[(partAResponses['ID#']
                                                                    == iD)
                                                                   & (partAResponses['Stim File']
                                                                      == 'A1.wav'),
                                                                   response].values)
        partAResponses.loc[(partAResponses['ID#'] == iD)
                           & (partAResponses['Stim File'].str.contains("A2")),
                           "d" + response] = (partAResponses.loc[(partAResponses['ID#']
                                                                  == iD)
                                                                 & (partAResponses['Stim File'].str.contains("A2"))].loc[:,
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

# change to highly annoyed data
partAData['dHighAnnoy'] = np.nan
for iD in partAData['ID#'].unique():
    
    # for participants with no HA in baseline, add HA occurences
    if partAData.loc[(partAData['ID#'] == iD) & (partAData['StimFile'] == 'A1.wav'), 'HighAnnoy'].iloc[0] == 0:
        partAData.loc[(partAData['ID#'] == iD)
                       & (partAData['StimFile'].str.contains("A1")),
                       'dHighAnnoy'] = partAData.loc[(partAData['ID#']
                                                      == iD)
                                                      & (partAData['StimFile'].str.contains("A1")),
                                                      'HighAnnoy']

    if partAData.loc[(partAData['ID#'] == iD) & (partAData['StimFile'] == 'A2.wav'), 'HighAnnoy'].iloc[0] == 0:
        partAData.loc[(partAData['ID#'] == iD)
                       & (partAData['StimFile'].str.contains("A2")),
                       'dHighAnnoy'] = partAData.loc[(partAData['ID#']
                                                      == iD)
                                                      & (partAData['StimFile'].str.contains("A2")),
                                                      'HighAnnoy']


# initialise DataFrames for loop over stimulus recordings
partA = pd.DataFrame(index=partAResponses['Recording'].unique())
partAstats = pd.DataFrame(index=partAResponses['Recording'].unique())

# loop over stimuli recording names, extract corresponding response data and
# tranpose data with recording as index
# apply basic normality tests and calculate aggregate statistics
for ii, file in enumerate(partAResponses['Recording'].unique()):
    if ii == 0:
        print("Processing results...\n")
    print(file.split('.')[0])
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
    
    partAdHighAnnoy = partAData.loc[partAData['Recording'] == file,
                                     ['ID#', 'dHighAnnoy']]
    columns = ["dHighAnnoy_" + str(ID) for ID in partAdHighAnnoy['ID#']]
    partAdHighAnnoy = pd.DataFrame(data=np.array(partAdHighAnnoy['dHighAnnoy']),
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
                                         'HighAnnoyProp'],
                                index=[file])

    dhighAnnoyAgg = pd.DataFrame(data=[[np.nansum(partAdHighAnnoy.values,
                                                  axis=1)[0],
                                        np.nanmean(partAdHighAnnoy.values,
                                                   axis=1)[0]]],
                                 columns=['dHighAnnoyTotal',
                                          'dHighAnnoyProp'],
                                 index=[file])

    noticeAgg = pd.DataFrame(data=[[np.sum(partANotice.values, axis=1)[0],
                                    np.mean(partANotice.values, axis=1)[0]]],
                             columns=['NoticedTotal', 'NoticedProp'],
                             index=[file])

    # add results to DataFrame
    if ii == 0:
        partA = partA.join([partAValence, partAArousal, partAAnnoy,
                            partAdValence, partAdArousal, partAdAnnoy,
                            partAHighAnnoy, partAdHighAnnoy,
                            partANotice, valenceAgg, arousalAgg,
                            annoyAgg, dvalenceAgg, darousalAgg,
                            dannoyAgg, highAnnoyAgg, dhighAnnoyAgg, noticeAgg])
        partAstats = partAstats.join([valenceSWtest, arousalSWtest,
                                      annoySWtest, dvalenceSWtest,
                                      darousalSWtest, dannoySWtest])
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
        partA.loc[file, partAdHighAnnoy.columns] = partAdHighAnnoy.loc[file]
        partA.loc[file, partANotice.columns] = partANotice.loc[file]
        partA.loc[file, valenceAgg.columns] = valenceAgg.loc[file]
        partA.loc[file, arousalAgg.columns] = arousalAgg.loc[file]
        partA.loc[file, annoyAgg.columns] = annoyAgg.loc[file]
        partA.loc[file, dvalenceAgg.columns] = dvalenceAgg.loc[file]
        partA.loc[file, darousalAgg.columns] = darousalAgg.loc[file]
        partA.loc[file, dannoyAgg.columns] = dannoyAgg.loc[file]
        partA.loc[file, highAnnoyAgg.columns] = highAnnoyAgg.loc[file]
        partA.loc[file, dhighAnnoyAgg.columns] = dhighAnnoyAgg.loc[file]
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

# change to highly annoyed data
partBData['dHighAnnoy'] = np.nan
for iD in partBData['ID#'].unique():

    if partBData.loc[(partBData['ID#'] == iD) & (partBData['StimFile'] == 'B2.wav'), 'HighAnnoy'].iloc[0] == 0:
        partBData.loc[(partBData['ID#'] == iD)
                       & (partBData['StimFile'].str.contains("B2")),
                       'dHighAnnoy'] = partBData.loc[(partBData['ID#']
                                                      == iD)
                                                      & (partBData['StimFile'].str.contains("B2")),
                                                      'HighAnnoy']

# initialise DataFrames for loop over stimulus recordings
partB = pd.DataFrame(index=partBResponses['Recording'].unique())

partBstats = pd.DataFrame(index=partBResponses['Recording'].unique())

# loop over stimuli recording names, extract corresponding response data and
# tranpose data with recording as index
# apply basic normality tests and calculate aggregate statistics
for ii, file in enumerate(partBResponses['Recording'].unique()):
    if ii == 0:
        print("Processing results...\n")
    print(file.split('.')[0])
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
    
    partBdHighAnnoy = partBData.loc[partBData['Recording'] == file,
                                     ['ID#', 'dHighAnnoy']]
    columns = ["dHighAnnoy_" + str(ID) for ID in partBdHighAnnoy['ID#']]
    partBdHighAnnoy = pd.DataFrame(data=np.array(partBdHighAnnoy['dHighAnnoy']),
                                   index=columns, columns=[file]).transpose()
    print(partBdHighAnnoy)
    
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
                                         'HighAnnoyProp'],
                                index=[file])

    dhighAnnoyAgg = pd.DataFrame(data=[[np.nansum(partBdHighAnnoy.values,
                                                  axis=1)[0],
                                        np.nanmean(partBdHighAnnoy.values,
                                                   axis=1)[0]]],
                                 columns=['dHighAnnoyTotal',
                                          'dHighAnnoyProp'],
                                 index=[file])

    # add results to DataFrame
    if ii == 0:
        partB = partB.join([partBValence, partBArousal, partBAnnoy,
                            partBdValence, partBdArousal, partBdAnnoy,
                            partBHighAnnoy, partBdHighAnnoy,
                            valenceAgg, arousalAgg, annoyAgg,
                            dvalenceAgg, darousalAgg, dannoyAgg, highAnnoyAgg,
                            dhighAnnoyAgg])
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
        partB.loc[file, partBdHighAnnoy.columns] = partBdHighAnnoy.loc[file]
        partB.loc[file, valenceAgg.columns] = valenceAgg.loc[file]
        partB.loc[file, arousalAgg.columns] = arousalAgg.loc[file]
        partB.loc[file, annoyAgg.columns] = annoyAgg.loc[file]
        partB.loc[file, dvalenceAgg.columns] = dvalenceAgg.loc[file]
        partB.loc[file, darousalAgg.columns] = darousalAgg.loc[file]
        partB.loc[file, dannoyAgg.columns] = dannoyAgg.loc[file]
        partB.loc[file, highAnnoyAgg.columns] = highAnnoyAgg.loc[file]
        partB.loc[file, dhighAnnoyAgg.columns] = dhighAnnoyAgg.loc[file]
    
allResponses = pd.concat([partB, partA], axis=0, join='outer')
allResponses.sort_index(inplace=True)
allResponses.set_index(allResponses.index.str.replace("_CALBIN_Pa.wav", ""),
                       inplace=True)

# merge Shapiro-Wilks test results
allShapWilksTest = pd.concat([partAstats, partBstats]).sort_index()
allShapWilksTest.set_index(allShapWilksTest.index.str.replace("_CALBIN_Pa.wav", ""),
                           inplace=True)

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

# ----------------------------------
# Prepare outputs for saving to file
# ----------------------------------

# separate 'by stimulus' output into test data and auxiliary data and save to
# file
dataByStimTest = dataByStim.loc[(dataByStim.index.str.find("B2") == 0)
                                | (dataByStim.index.str.find("A1") == 0)
                                | (dataByStim.index.str.find("A2") == 0), :]
dataByStimTestA = dataByStimTest.loc[(dataByStimTest.index.str.find("A1") == 0)
                                     | (dataByStimTest.index.str.find("A2") == 0),
                                     :].dropna(axis=1, how='all')
dataByStimTestB = dataByStimTest.loc[dataByStimTest.index.str.find("B2") == 0,
                                     :].dropna(axis=1, how='all')

dataByStimAux = dataByStim.loc[~((dataByStim.index.str.find("B2") == 0)
                                 | (dataByStim.index.str.find("A1") == 0)
                                 | (dataByStim.index.str.find("A2") == 0)),
                               :].copy().dropna(axis=1, how='all')

# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

outFilePath = QFileDialog.getExistingDirectory(caption="Choose output folder to save processed files in '03 Experiment\Experiment 1\Analysis\PostProcess'")


dataByStim.to_csv(os.path.join(outFilePath,
                               "refmap_listest1_alldata_ByStim.csv"))

dataByStimTest.to_csv(os.path.join(outFilePath,
                                   "refmap_listest1_testdata_ByStim.csv"))

dataByStimTestA.to_csv(os.path.join(outFilePath,
                                    "refmap_listest1_testdataA_ByStim.csv"))

dataByStimTestB.to_csv(os.path.join(outFilePath,
                                    "refmap_listest1_testdataB_ByStim.csv"))

dataByStimAux.to_csv(os.path.join(outFilePath,
                                  "refmap_listest1_auxdata.csv"))

# save statistical tests
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

statsoutFilePath = QFileDialog.getExistingDirectory(caption="Choose output folder to save statistical analysis files in '03 Experiment\Experiment 1\Analysis\Python'")
allShapWilksTest.to_csv(os.path.join(statsoutFilePath,
                                     "refmap_listest1_alltestShapWilks.csv"))

# merge response and stimuli data into 'by participant' test datasets, and
# save to file

partADataBySubj = pd.merge(left=partAData,
                           right=dataByStimTestA.loc[:, :indicesDiffPsycho[-1]],
                           how='outer', left_on='Recording', right_on='CALBINRecFiles')
partADataBySubj.sort_values(by='ID#', axis=0, inplace=True)
partADataBySubj = pd.merge(left=partADataBySubj,
                           right=prePostTestResponses, how='left',
                           left_on='ID#', right_on='ID#')
partADataBySubj.drop(columns=['Part', 'Recording'], inplace=True)
partADataBySubj.insert(loc=0, column='SessionPart',
                       value=partADataBySubj.pop('SessionPart'))

partBDataBySubj = pd.merge(left=partBData,
                           right=dataByStimTestB.loc[:, :indicesDiffPsycho[-1]],
                           how='left', left_on='Recording', right_on='CALBINRecFiles')
partBDataBySubj.sort_values(by='ID#', axis=0, inplace=True)
partBDataBySubj = pd.merge(left=partBDataBySubj,
                           right=prePostTestResponses, how='outer',
                           left_on='ID#', right_on='ID#')
partBDataBySubj.drop(columns=['Part', 'Recording'], inplace=True)
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
omitColumns = omitColumns + (["dHighAnnoy_" + str(partID) for partID in omitParticipants])
dataByStimTestANotice = dataByStimTestA.drop(labels=omitColumns, axis=1)
dataByStimTestANotice.drop(labels=['ArousalMean', 'ArousalMedian',
                                   'ValenceMean', 'ValenceMedian',
                                   'AnnoyMean', 'AnnoyMedian',
                                   'dArousalMean', 'dArousalMedian',
                                   'dValenceMean', 'dValenceMedian',
                                   'dAnnoyMean', 'dAnnoyMedian',
                                   'HighAnnoyTotal', 'HighAnnoyProp',
                                   'dHighAnnoyTotal', 'dHighAnnoyProp',
                                   'NoticedTotal', 'NoticedProp'], axis=1,
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
keepColumns = [label.replace("HighAnnoy_", "dHighAnnoy_") for label in keepColumns]
dataByStimTestANotice['dHighAnnoyTotalFilt'] = dataByStimTestANotice[keepColumns].sum(axis=1)
dataByStimTestANotice['dHighAnnoyPropFilt'] = dataByStimTestANotice['dHighAnnoyTotalFilt']/len(keepParticipants)

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
                                                 columns='StimFile',
                                                 values='Annoyance')
partAValenceDataBySubjWide = partADataBySubj.pivot(index='ID#',
                                                   columns='StimFile',
                                                   values='Valence')
partAArousalDataBySubjWide = partADataBySubj.pivot(index='ID#',
                                                   columns='StimFile',
                                                   values='Arousal')
partAdAnnoyDataBySubjWide = partADataBySubj.pivot(index='ID#',
                                                  columns='StimFile',
                                                  values='dAnnoyance')
partAdValenceDataBySubjWide = partADataBySubj.pivot(index='ID#',
                                                    columns='StimFile',
                                                    values='dValence')
partAdArousalDataBySubjWide = partADataBySubj.pivot(index='ID#',
                                                    columns='StimFile',
                                                    values='dArousal')
partANoticeDataBySubjWide = partADataBySubj.pivot(index='ID#',
                                                  columns='StimFile',
                                                  values='UAS_noticed')
partANoticeFiltDataBySubjWide = partADataBySubjNotice.pivot(index='ID#',
                                                            columns='StimFile',
                                                            values='UAS_noticed')

# merge all together for multivariate analysis
partADataBySubjWide = partADataBySubj.pivot(index='ID#',
                                            columns='StimFile',
                                            values='Annoyance')
partADataBySubjWide.columns = partADataBySubjWide.columns.str.replace(".wav", "_annoy")
partADataBySubjWide = partADataBySubjWide.merge(partAValenceDataBySubjWide,
                                                on='ID#')
partADataBySubjWide.columns = partADataBySubjWide.columns.str.replace(".wav", "_valence")
partADataBySubjWide = partADataBySubjWide.merge(partAArousalDataBySubjWide,
                                                on='ID#')
partADataBySubjWide.columns = partADataBySubjWide.columns.str.replace(".wav", "_arousal")
partADataBySubjWide = partADataBySubjWide.merge(partANoticeDataBySubjWide,
                                                on='ID#')
partADataBySubjWide.columns = partADataBySubjWide.columns.str.replace(".wav", "_notice")

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
                                                 columns='StimFile',
                                                 values='Annoyance')
partBValenceDataBySubjWide = partBDataBySubj.pivot(index='ID#',
                                                   columns='StimFile',
                                                   values='Valence')
partBArousalDataBySubjWide = partBDataBySubj.pivot(index='ID#',
                                                   columns='StimFile',
                                                   values='Arousal')
partBdAnnoyDataBySubjWide = partBDataBySubj.pivot(index='ID#',
                                                  columns='StimFile',
                                                  values='dAnnoyance')
partBdValenceDataBySubjWide = partBDataBySubj.pivot(index='ID#',
                                                    columns='StimFile',
                                                    values='dValence')
partBdArousalDataBySubjWide = partBDataBySubj.pivot(index='ID#',
                                                    columns='StimFile',
                                                    values='dArousal')

# merge all together for multivariate analysis
partBDataBySubjWide = partBDataBySubj.pivot(index='ID#',
                                            columns='StimFile',
                                            values='Annoyance')
partBDataBySubjWide.columns = partBDataBySubjWide.columns.str.replace(".wav", "_annoy")
partBDataBySubjWide = partBDataBySubjWide.merge(partBValenceDataBySubjWide,
                                                on='ID#')
partBDataBySubjWide.columns = partBDataBySubjWide.columns.str.replace(".wav", "_valence")
partBDataBySubjWide = partBDataBySubjWide.merge(partBArousalDataBySubjWide,
                                                on='ID#')
partBDataBySubjWide.columns = partBDataBySubjWide.columns.str.replace(".wav", "_arousal")

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