# -*- coding: utf-8 -*-

# script


# --------
# %% Setup
# --------

# import statements
import sys
import os
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import QFileDialog, QApplication
import librosa
from refmap_psychoacoustics.dsp import filterFuncs
from refmap_psychoacoustics.metrics import psych_annoy
from scipy import stats, io
from warnings import simplefilter

# suppress pandas performance warnings
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

# enable copy-on-write mode for Pandas (will be default from Pandas 3.0)
pd.options.mode.copy_on_write = True

# ------------------------------------------------------------
# %% Initialise data frame for metric calculations and results
# ------------------------------------------------------------

# skip signal start and end for time-aggregation
start_skipT = 1.5
end_skipT = 1.5

# open wav file selection dialog and assign filepaths to list
# PROJECT NOTE: the calibrated HATS files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%202/Stimuli/Calibrated_recordings/RecordHATS?csf=1&web=1&e=WNwAvH
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance()

fileExts = "*.wav"
filelist = list(QFileDialog.getOpenFileNames(caption="Open recording files in '03 Experiment\Experiment 2\Calibration\Stimuli\Calibrated_recordings\RecordHATS'",
                                             filter=fileExts))[0]
filelist.sort()
filenames = [filepath.split('/')[-1] for filepath in filelist]
stemNames = [filename.replace("_HATS_Pa.wav", "") for filename in filenames]
sqmNames = [filename.replace("HATS", "MA220") for filename in filenames
            if filename.find("G57") == -1]


dataByStim = pd.DataFrame(index=stemNames, dtype=float)

dataByStim['HATSRecFiles'] = filenames
dataByStim['MA220MicRecFiles'] = sqmNames

# ----------------------------------------
# %% Add categorical variables for stimuli
# ----------------------------------------

# NOTE: the order of this section matters, as some categorical variable
# definitions depend on Boolean logic derived from other categorical variables
# added

# give each test stimulus (only those prefixed by an ambient env label) a unique ID number that is different to all other stimuli
testStimMask = ((dataByStim.index.str.find("Street") != -1)
                | (dataByStim.index.str.find("Park") != -1))
testStimIndices = dataByStim.index[testStimMask]
for ii, file in enumerate(testStimIndices):
    dataByStim.loc[file, 'StimID'] = ii + 1
dataByStim['StimID'] = dataByStim['StimID'].astype('Int64')

# ambient env category
dataByStim.loc[dataByStim.index.str.find("Background") == -1, 'AmbientEnv'] = "UAS only"
dataByStim.loc[dataByStim.index.str.find("BusyStreet6") != -1, 'AmbientEnv'] = "Highway"
dataByStim.loc[dataByStim.index.str.find("BusyStreet6Quiet") != -1, 'AmbientEnv'] = "Highway (low)"
dataByStim.loc[dataByStim.index.str.find("BusyStreet8") != -1, 'AmbientEnv'] = "Streetside square"
dataByStim.loc[dataByStim.index.str.find("Park3") != -1, 'AmbientEnv'] = "Park"
dataByStim.loc[dataByStim.index.str.find("Park3Loud") != -1, 'AmbientEnv'] = "Park (high)"
dataByStim.loc[dataByStim.index.str.find("QuietStreet7") != -1, 'AmbientEnv'] = "Residential"

# UAS operation
dataByStim.loc[dataByStim.index.str.find("Overflight") != -1, 'UASOperation'] = "Overflight"
dataByStim.loc[dataByStim.index.str.find("Delivery") != -1, 'UASOperation'] = "Delivery"
dataByStim.loc[(dataByStim.index.str.find("Baseline") != -1),
               'UASOperation'] = "Baseline"

# UAS event quantity
dataByStim.loc[dataByStim.index.str.find("Baseline") != -1, 'UASEvents'] = 0
dataByStim.loc[dataByStim.index.str.find("_1") != -1, 'UASEvents'] = 1
dataByStim.loc[dataByStim.index.str.find("_2") != -1, 'UASEvents'] = 2
dataByStim.loc[dataByStim.index.str.find("_3") != -1, 'UASEvents'] = 3
dataByStim['UASEventPerMin'] = dataByStim['UASEvents']/(30/60)

# UAS type
dataByStim.loc[dataByStim.index.str.find("H520") != -1, 'UASType'] = "H520"
dataByStim.loc[dataByStim.index.str.find("T150") != -1, 'UASType'] = "T150"
dataByStim.loc[dataByStim.index.str.find("Baseline") != -1, 'UASType'] = "Baseline"

# UAS starting hemisphere (side)
dataByStim.loc[dataByStim.index.str.find("left") != -1, 'UASStart'] = "Left"
dataByStim.loc[dataByStim.index.str.find("right") != -1, 'UASStart'] = "Right"

# --------------------------------------------------
# %% Acoustic, PNL & Detection metrics single values
# --------------------------------------------------

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
    signalA = filterFuncs.A_weight_T(signal, sampleRatein)
    signalmagAF = filterFuncs.time_weight(signalA, sampleRatein, tau=0.125)
    signalmagAS = filterFuncs.time_weight(signalA, sampleRatein, tau=1)

    # calculate weighted dB time series
    signaldBAF = 20*np.log10(signalmagAF[start_skips:-end_skips]/2e-5)
    signaldBAS = 20*np.log10(signalmagAS[start_skips:-end_skips]/2e-5)

    # calculate percentile metrics
    signalLAFmax = signaldBAF.max(axis=0)
    signalLAF5 = np.percentile(signaldBAF, q=95, axis=0, method='median_unbiased')
    signalLAF10 = np.percentile(signaldBAF, q=90, axis=0, method='median_unbiased')
    signalLAF25 = np.percentile(signaldBAF, q=75, axis=0, method='median_unbiased')
    signalLAF50 = np.percentile(signaldBAF, q=50, axis=0, method='median_unbiased')
    signalLAF75 = np.percentile(signaldBAF, q=25, axis=0, method='median_unbiased')
    signalLAF90 = np.percentile(signaldBAF, q=10, axis=0, method='median_unbiased')
    signalLAF95 = np.percentile(signaldBAF, q=5, axis=0, method='median_unbiased')
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

    if dataByStim.loc[stemNames[ii], 'AmbientEnv'] != "UAS only":

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


# end of for loop over HATS signal wav files


# -----------------------------------
# %% PNL and detection metrics import
# -----------------------------------

indicesPNL = ["PNLmaxMaxLR", "PNLTmaxMaxLR", "EPNLMaxLR"]
indicesDetect = ["Detect0p5dBADiscMaxLR", "Detect0p5dBMaxMaxLR",
                 "Detect0p5dBIntMaxLR", "Detect0p1dBADiscMaxLR",
                 "Detect0p1dBMaxMaxLR", "Detect0p1dBIntMaxLR",
                 "Detect0p5dBMaxEx50MaxLR", "Detect0p5dBIntEx50MaxLR",
                 "Detect0p1dBMaxEx50MaxLR", "Detect0p1dBIntEx50MaxLR"]

# import PNL and detection results
fileExts = "*.csv"
filelist = list(QFileDialog.getOpenFileNames(caption="Open results files in '03 Experiment\Experiment 2\Analysis\MATLAB\HATS\csv'",
                                             filter=fileExts))[0]
filelist.sort()
# assumes read files in alphabetical order with NASA_ first and SQAT_ second
DetectResults = pd.read_csv(filelist[0], header=0, index_col=0)
PNLResults = pd.read_csv(filelist[1], header=0, index_col=0)
PNLDetectResults = PNLResults.merge(DetectResults, left_index=True, right_index=True)
PNLDetectResults.set_index(PNLDetectResults.index.str.replace(pat="_HATS_Pa.wav", repl=""), inplace=True)

# join results to data
dataByStim = dataByStim.join(PNLDetectResults, how='left')

# ---------------------------------------------------
# %% Psychoacoustic metrics single values calculation
# ---------------------------------------------------

# output variables
indicesPsycho = ["LoudECMAPowAvg",
                 "LoudISO105Ex",
                 "LoudISO1PowAvg",
                 "LoudISO3PowAvg",
                 "LoudISO305Ex",
                 "LoudQZ5321PowAvg",
                 "LoudQZ532105Ex",
                 "LoudQZ5323PowAvg",
                 "LoudQZ532305Ex",
                 "LoudQZ4182PowAvg",
                 "LoudQZ418205Ex",
                 "TonalECMAAvg",
                 "TonalECMA05Ex",
                 "TonalSHMIntAvg",
                 "TonalSHMInt05Ex",
                 "TonalAwSHMAvg",
                 "TonalAwSHM05Ex",
                 "TonalAwSHMIntAvg",
                 "TonalAwSHMInt05Ex",
                 "TonLdECMAPowAvg",
                 "TonLdECMA05Ex",
                 "TonShpAurSHMPowAvg",
                 "TonShpAurSHM05Ex",
                 "TonalAurAvg",
                 "TonalAur10Ex",
                 "TonalAur05Ex",
                 "RoughECMA10Ex",
                 "RoughECMA05Ex",
                 "RoughECMAAvg",
                 "RoughFZ10Ex",
                 "RoughFZ05Ex",
                 "RoughDW10Ex",
                 "RoughDW05Ex",
                 "FluctOldSHM10Ex",
                 "FluctOldSHM05Ex",
                 "FluctOldSHMAvg",
                 "FluctECMA10Ex",
                 "FluctECMA05Ex",
                 "FluctFZ10Ex",
                 "FluctFZ05Ex",
                 "FluctOV10Ex",
                 "FluctOV05Ex",
                 "SharpAurSHMPowAvg",
                 "SharpAurSHM05Ex",
                 "SharpAurISO3PowAvg",
                 "SharpAurISO305Ex",
                 "SharpAurISO1PowAvg",
                 "SharpAurISO105Ex",
                 "SharpAurISO1Med",
                 "SharpDINPowAvg",
                 "SharpDIN05Ex",
                 "SharpvBISO1PowAvg",
                 "SharpvBISO105Ex",
                 "SharpAurQZ5321PowAvg",
                 "SharpAurQZ532105Ex",
                 "SharpAurQZ5323PowAvg",
                 "SharpAurQZ532305Ex",
                 "SharpAurQZ4182PowAvg",
                 "SharpAurQZ418205Ex",
                 "ImpulsSHMPowAvg",
                 "ImpulsSHMAvg",
                 "ImpulsSHM05Ex",
                 "ImpulsLoudWZAvg",
                 "ImpulsLoudWZ05Ex",
                 "ImpulsLoudWZPowAvg",
                 "ImpulsLoudWECMAAvg",
                 "ImpulsLoudWECMA05Ex",
                 "ImpulsLoudWECMAPowAvg",
                 "PsychAnnoyWidmann",
                 "PsychAnnoyMore",
                 "PsychAnnoyDi",
                 "PsychAnnoyTorija",
                 "PsychAnnoyWillemsen",
                 "PsychAnnoyBoucher",
                 "PartLoudSHMPowAvg",
                 "PartTonLdSHMPowAvg",
                 "PartTonShpAurSHMPowAvg",
                 "PartTonShpAurSHM05Ex",
                 "dTonalECMAAvg",
                 "dTonalECMA05Ex",
                 "dTonalSHMIntAvg",
                 "dTonalSHMInt05Ex",
                 "dTonalAwSHMAvg",
                 "dTonalAwSHM05Ex",
                 "dTonalAwSHMIntAvg",
                 "dTonalAwSHMInt05Ex",
                 "dTonLdECMAPowAvg",
                 "dTonLdECMA05Ex",
                 "dTonShpAurSHMPowAvg",
                 "dTonShpAurSHM05Ex",
                 "dRoughECMA10Ex",
                 "dRoughECMA05Ex",
                 "dRoughFZ10Ex",
                 "dRoughFZ05Ex",
                 "dFluctECMA10Ex",
                 "dFluctECMA05Ex",
                 "dFluctOV10Ex",
                 "dFluctOV05Ex",
                 "dSharpAurSHMPowAvg",
                 "dSharpAurSHM05Ex",
                 "dSharpAurISO3PowAvg",
                 "dSharpAurISO305Ex",
                 "dImpulsSHMPowAvg",
                 "dImpulsSHMAvg",
                 "dImpulsSHM05Ex",
                 "dImpulsLoudWZAvg",
                 "dImpulsLoudWZ05Ex",
                 "dImpulsLoudWZPowAvg",
                 "dImpulsLoudWECMAAvg",
                 "dImpulsLoudWECMA05Ex",
                 "dImpulsLoudWECMAPowAvg",
                 "dPsychAnnoyWidmann",
                 "dPsychAnnoyMore",
                 "dPsychAnnoyDi",
                 "dPsychAnnoyTorija",
                 "dPsychAnnoyWillemsen",
                 "dPsychAnnoyBoucher"]

dataByStim = pd.concat([dataByStim, pd.DataFrame(index=dataByStim.index,
                                                 columns=indicesPsycho,
                                                 dtype=float)], axis=1)


# output SQM sample rates
sampleRateLoudECMA = 187.5
sampleRateLoudISO1 = 500
sampleRateLoudISO3 = 1e3
sampleRateTonalECMA = 187.5
sampleRateTonalAures = 12.5
sampleRateModECMA = 50
sampleRateRoughFZ = 2000
sampleRateRoughDW = 10
sampleRateFluctOldSHM = 229390681/4e6
sampleRateFluctOV = 5
sampleRateImpulsSHM = 60000/55
sampleRateQZ = 10
sampleRateLoudISO1HighRes = 2000

# critical band differences (used for integration of specific values)
bandDiff0p5 = 0.5  # used for all except Sottek Hearing Model impulsiveness
bandDiffSHMImp = 1.0  # Sottek Hearing Model impulsiveness band differences

# time value for moving averaging of SQMs for difference calculations
windowT = 0.05

# From MATLAB calculations
# ------------------------

# open xlsx file selection dialog and assign filepaths to list
# PROJECT NOTE: the calculated files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%202/Analysis/MATLAB/MA220Mic?csf=1&web=1&e=23deLX
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

fileExtsXLSX = "*.xlsx"
filelistXLSX = list(QFileDialog.getOpenFileNames(filter=fileExtsXLSX,
                                                 caption=r"Select MATLAB output files in '03 Experiment\Experiment 2\Analysis\MATLAB\MA220Mic\xlsx'"))[0]
filelistXLSX.sort()

fileExtsMAT = "*.mat"
filelistMAT = list(QFileDialog.getOpenFileNames(filter=fileExtsMAT,
                                                caption=r"Select MATLAB output files in '03 Experiment\Experiment 2\Analysis\MATLAB\MA220Mic\mat'"))[0]

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
        if filenames[ii].find("FluctFZ") != -1:
            filedata = io.loadmat(file)['FluctFZSpecTDep']

            # Calculate overall Fastl & Zwicker fluctuation strength from specific
            # time-dependent fluctuation strength
            specFluctFZ = pd.DataFrame(filedata[1:, 1:],
                                       index=filedata[1:, 0],
                                       columns=filedata[0, 1:])

            fluctFZTDep = bandDiff0p5*specFluctFZ.sum(axis=1).to_frame(name="Monaural")

            # mask for start/end skip and 0 values
            fluctFZTDepMask = fluctFZTDep.loc[(fluctFZTDep.index.values
                                               > start_skipT).transpose()
                                              & (fluctFZTDep.index.values
                                                 < fluctFZTDep.index.values.max()
                                                 - end_skipT).transpose()]

            # overall 5% exceeded fluctuation strength
            fluctFZ05Ex = np.percentile(fluctFZTDepMask, q=95, axis=0, 
                                        method='median_unbiased')

            # overall 10% exceeded fluctuation strength
            fluctFZ10Ex = np.percentile(fluctFZTDepMask, q=90, axis=0,
                                        method='median_unbiased')
            
            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'FluctFZ10Ex'] = fluctFZ10Ex
            dataByStim.loc[renderNames[ii], 'FluctFZ05Ex'] = fluctFZ05Ex

        elif filenames[ii].find("RoughFZ") != -1:
            filedata = io.loadmat(file)['RoughFZSpecTDep']

            # Calculate overall Fastl & Zwicker roughness from specific
            # time-dependent roughness
            specRoughFZ = pd.DataFrame(filedata[1:, 1:],
                                        index=filedata[1:, 0],
                                        columns=filedata[0, 1:])

            roughFZTDep = bandDiff0p5*specRoughFZ.sum(axis=1).to_frame(name="Monaural")

            # mask for start/end skip and 0 values
            roughFZTDepMask = roughFZTDep.loc[(roughFZTDep.index.values
                                               > start_skipT).transpose()
                                              & (roughFZTDep.index.values
                                                 < roughFZTDep.index.values.max()
                                                 - end_skipT).transpose()]

            # overall 5% exceeded roughness
            roughFZ05Ex = np.percentile(roughFZTDepMask, q=95, axis=0, 
                                        method='median_unbiased')

            # overall 10% exceeded roughness
            roughFZ10Ex = np.percentile(roughFZTDepMask, q=90, axis=0,
                                        method='median_unbiased')

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'RoughFZ10Ex'] = roughFZ10Ex
            dataByStim.loc[renderNames[ii], 'RoughFZ05Ex'] = roughFZ05Ex

        elif filenames[ii].find("LoudECMA") != -1:
            filedata = io.loadmat(file)['LoudECMASpecTDep']

            # Calculate ECMA-418-2:2025 Sottek Hearing Model overall loudness from
            # specific time-dependent loudness
            specLoudECMA = pd.DataFrame(filedata[1:, 1:],
                                         index=filedata[1:, 0],
                                         columns=filedata[0, 1:])

            # time-dependent loudness (ECMA-418-2:2025 Equation 116)
            loudECMATDep = specLoudECMA.sum(axis=1)*bandDiff0p5

            # mask for start/end skip
            loudECMATDepMask = loudECMATDep.loc[(loudECMATDep.index.values
                                                 > start_skipT)
                                                & (loudECMATDep.index.values
                                                   < loudECMATDep.index.values.max()
                                                   - end_skipT)]

            # overall (power-averaged) loudness (ECMA-418-2:2022 Equation 117)
            loudECMAPowAvg = loudECMATDepMask.pow(1/np.log10(2)).mean()**np.log10(2)

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'LoudECMAPowAvg'] = loudECMAPowAvg

        elif filenames[ii].find("TonLdECMA") != -1:
            filedata = io.loadmat(file)['TonLdECMASpecTDep']

            # Calculate ECMA-418-2:2025 Sottek Hearing Model overall tonal loudness from
            # specific time-dependent tonal loudness
            specTonLdECMA = pd.DataFrame(filedata[1:, 1:],
                                          index=filedata[1:, 0],
                                          columns=filedata[0, 1:])

            # time-dependent tonal loudness (ECMA-418-2:2025 Equation 116)
            tonLdECMATDep = specTonLdECMA.sum(axis=1)*bandDiff0p5

            # mask for start/end skip
            tonLdECMATDepMask = tonLdECMATDep.loc[(tonLdECMATDep.index.values
                                                   > start_skipT)
                                                  & (tonLdECMATDep.index.values
                                                     < tonLdECMATDep.index.values.max()
                                                     - end_skipT)]

            # overall (power-averaged) tonal loudness (ECMA-418-2:2025
            # Equation 117)
            tonLdECMAPowAvg = tonLdECMATDepMask.pow(1/np.log10(2)).mean()**np.log10(2)

            # 5% exceeded tonal loudness
            tonLdECMA05Ex = np.percentile(tonLdECMATDepMask, q=95, axis=0,
                                          method='median_unbiased')
            
            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'TonLdECMAPowAvg'] = tonLdECMAPowAvg
            dataByStim.loc[renderNames[ii], 'TonLdECMA05Ex'] = tonLdECMA05Ex
        
        elif filenames[ii].find("TonalECMA") != -1:
            filedata = io.loadmat(file)['TonalECMASpecTDep']

            # Calculate ECMA-418-2:2025 Sottek Hearing Model overall tonality from
            # specific tonality
            specTonalECMA = pd.DataFrame(filedata[1:, 1:],
                                         index=filedata[1:, 0],
                                         columns=filedata[0, 1:])

            # time-dependent tonality (max, not integration)
            tonalECMATDep = specTonalECMA.max(axis=1)

            # time-dependent tonality (integrated, with recalibration adjustment
            # to match 40 dB 1 kHz sine to 1 tu)
            tonalSHMIntTDep = specTonalECMA.sum(axis=1)*bandDiff0p5*0.348088948583815

            # mask for start/end skip and values <= 0.02
            tonalECMATDepMask = tonalECMATDep.loc[(tonalECMATDep.index.values
                                                   > start_skipT)
                                                  & (tonalECMATDep.index.values
                                                     < tonalECMATDep.index.values.max()
                                                     - end_skipT)
                                                  & (tonalECMATDep.values
                                                     > 0.02)]

            # mask for start/end skip and values <= 0.02
            # NOTE: uses mask from ECMA tonality
            tonalSHMIntTDepMask = tonalSHMIntTDep.loc[(tonalSHMIntTDep.index.values
                                                         > start_skipT)
                                                        & (tonalSHMIntTDep.index.values
                                                           < tonalSHMIntTDep.index.values.max()
                                                           - end_skipT)
                                                        & (tonalECMATDep.values
                                                           > 0.02)]  # see NOTE above

            # time-averaged tonality (omitting T<=0.02)
            tonalECMAAvg = tonalECMATDepMask.mean(axis=0)

            # 5% exceeded tonality (omitting T<=0.02)
            tonalECMA05Ex = np.percentile(tonalECMATDepMask, q=95, axis=0,
                                          method='median_unbiased')

            # time-averaged integrated tonality (omitting T<=0.02)
            # NOTE: uses mask from ECMA tonality
            tonalSHMIntAvg = tonalSHMIntTDepMask.mean(axis=0)

            # max l/r 5% exceeded integrated tonality (omitting T<=0.02)
            tonalSHMInt05Ex = np.percentile(tonalSHMIntTDepMask, q=95, axis=0,
                                            method='median_unbiased')

            # add results to output DataFrame            
            dataByStim.loc[renderNames[ii], 'TonalECMAAvg'] = tonalECMAAvg
            dataByStim.loc[renderNames[ii], 'TonalECMA05Ex'] = tonalECMA05Ex
            dataByStim.loc[renderNames[ii], 'TonalSHMIntAvg'] = tonalSHMIntAvg
            dataByStim.loc[renderNames[ii], 'TonalSHMInt05Ex'] = tonalSHMInt05Ex

            # Calculate tonal annoyance-weighted tonality from specific
            # time-dependent tonality
            band_centre_freqs = (81.9289/0.1618)*np.sinh(0.1618*np.arange(0.5, 27, 0.5))
            annoyWeight = np.ones(band_centre_freqs.shape)
            annoyWeight[band_centre_freqs
                        >= 1e3] = 2.3*(np.log10(band_centre_freqs[band_centre_freqs
                                                                  >= 1e3])
                                                                  - 3) + 1
            
            specTonalAwSHM = specTonalECMA.multiply(annoyWeight, axis=1)

            # time-dependent tonal annoyance-weighted tonality
            # (max, not integration)
            tonalAwSHMTDep = specTonalAwSHM.max(axis=1)

            # time-dependent tonal annoyance-weighted tonality
            # (integrated, with recalibration adjustment)
            tonalAwSHMIntTDep = specTonalAwSHM.sum(axis=1)*bandDiff0p5*0.348088948583815

            # mask for start/end skip and values <= 0.02
            tonalAwSHMTDepMask = tonalAwSHMTDep.loc[(tonalAwSHMTDep.index.values
                                                     > start_skipT)
                                                    & (tonalAwSHMTDep.index.values
                                                       < tonalAwSHMTDep.index.values.max()
                                                         - end_skipT)
                                                      & (tonalAwSHMTDep.values
                                                         > 0.02)]

            # mask for start/end skip and values <= 0.02
            # NOTE: uses mask from ECMA tonal annoyance-weighted tonality
            tonalAwSHMIntTDepMask = tonalAwSHMIntTDep.loc[(tonalAwSHMIntTDep.index.values
                                                           > start_skipT)
                                                          & (tonalAwSHMIntTDep.index.values
                                                             < tonalAwSHMIntTDep.index.values.max()
                                                             - end_skipT)
                                                          & (tonalAwSHMTDep.values
                                                             > 0.02)]  # see NOTE above

            # time-averaged tonal annoyance-weighted tonality (omitting T<=0.02)
            tonalAwSHMAvg = tonalAwSHMTDepMask.mean(axis=0)

            # 5% exceeded tonal annoyance-weighted tonality (omitting T<=0.02)
            tonalAwSHM05Ex = np.percentile(tonalAwSHMTDepMask, q=95, axis=0,
                                           method='median_unbiased')

            # time-averaged integrated tonal annoyance-weighted tonality (omitting T<=0.02)
            # NOTE: uses mask from ECMA tonal annoyance-weighted tonality
            tonalAwSHMIntAvg = tonalAwSHMIntTDepMask.mean(axis=0)

            # 5% exceeded integrated tonal annoyance-weighted tonality (omitting T<=0.02)
            tonalAwSHMInt05Ex = np.percentile(tonalAwSHMIntTDepMask, q=95, axis=0,
                                              method='median_unbiased')

            # add results to output DataFrame            
            dataByStim.loc[renderNames[ii], 'TonalAwSHMAvg'] = tonalAwSHMAvg
            dataByStim.loc[renderNames[ii], 'TonalAwSHM05Ex'] = tonalAwSHM05Ex
            dataByStim.loc[renderNames[ii], 'TonalAwSHMIntAvg'] = tonalAwSHMIntAvg
            dataByStim.loc[renderNames[ii], 'TonalAwSHMInt05Ex'] = tonalAwSHMInt05Ex

        elif filenames[ii].find("RoughECMA") != -1:
            filedata = io.loadmat(file)['RoughECMASpecTDep']

            # Calculate ECMA-418-2:2022 Sottek Hearing Model overall roughness from
            # specific time-dependent roughness
            specRoughECMA = pd.DataFrame(filedata[1:, 1:],
                                         index=filedata[1:, 0],
                                         columns=filedata[0, 1:])

            # time-dependent roughness (ECMA-418-2:2022 Equation 110)
            roughECMATDep = specRoughECMA.sum(axis=1)*bandDiff0p5

            # mask for start/end skip
            roughECMATDepMask = roughECMATDep.loc[(roughECMATDep.index.values
                                                   > start_skipT)
                                                  & (roughECMATDep.index.values
                                                     < roughECMATDep.index.values.max()
                                                     - end_skipT)]

            # overall (90th percentile = 10% exceeded) roughness
            roughECMA10Ex = np.percentile(roughECMATDepMask, q=90, axis=0,
                                          method='median_unbiased')
            # overall (95th percentile = 5% exceeded) roughness
            roughECMA05Ex = np.percentile(roughECMATDepMask, q=95, axis=0,
                                          method='median_unbiased')

            # overall time-averaged roughness
            roughECMAAvg = roughECMATDepMask.mean()

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'RoughECMA10Ex'] = roughECMA10Ex
            dataByStim.loc[renderNames[ii], 'RoughECMA05Ex'] = roughECMA05Ex
            dataByStim.loc[renderNames[ii], 'RoughECMAAvg'] = roughECMAAvg

        elif filenames[ii].find("FluctOV") != -1:       
            filedata = io.loadmat(file)['FluctOVSpecTDep']   
  
            # Calculate Osses Vecchi et al overall fluctuation strength from
            # specific time-dependent fluctuation strength
            specFluctOV = pd.DataFrame(filedata[1:, 1:],
                                       index=filedata[1:, 0],
                                       columns=filedata[0, 1:])

            # time-dependent fluctuation strength
            fluctOVTDep = specFluctOV.sum(axis=1)*bandDiff0p5

            # mask for start/end skip
            fluctOVTDepMask = fluctOVTDep.loc[(fluctOVTDep.index.values
                                               > start_skipT).transpose()
                                              & (fluctOVTDep.index.values
                                                 < fluctOVTDep.index.values.max()
                                                 - end_skipT).transpose()]

            # overall (90th percentile = 10% exceeded) fluctuation strength
            fluctOV10Ex = np.percentile(fluctOVTDepMask, q=90, axis=0,
                                        method='median_unbiased')

            # overall (95th percentile = 5% exceeded) fluctuation strength
            fluctOV05Ex = np.percentile(fluctOVTDepMask, q=95, axis=0,
                                        method='median_unbiased')

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'FluctOV10Ex'] = fluctOV10Ex
            dataByStim.loc[renderNames[ii], 'FluctOV05Ex'] = fluctOV05Ex
        
        elif filenames[ii].find("ImpulsLoudWZ") != -1:
            filedata = io.loadmat(file)['ImpulsLoudWZTDep']
            
            # Calculate Willemen et al overall impulsive loudness from 
            # time-dependent values
            impulsLoudWZTDep = pd.DataFrame(filedata[1:, 1:],
                                            index=filedata[1:, 0],
                                            columns=["Monaural"])
            # mask for start/end skip
            impulsLoudWZTDepMask = impulsLoudWZTDep.loc[(impulsLoudWZTDep.index.values
                                                         > start_skipT)
                                                        & (impulsLoudWZTDep.index.values
                                                           < impulsLoudWZTDep.index.values.max()
                                                           - end_skipT), :]
            
            # overall (averaged) impulsive loudness
            impulsLoudWZAvg = impulsLoudWZTDepMask.mean(axis=0).iloc[0]

            # overall (95th percentile = 5% exceeded) impulsive loudness
            impulsLoudWZ05Ex = np.percentile(impulsLoudWZTDepMask, q=95, axis=0,
                                             method='median_unbiased')

            # overall (power-averaged) impulsive loudness
            impulsLoudWZPowAvg = (impulsLoudWZTDepMask.pow(1/np.log10(2)).mean(axis=0)**np.log10(2)).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'ImpulsLoudWZAvg'] = impulsLoudWZAvg
            dataByStim.loc[renderNames[ii], 'ImpulsLoudWZ05Ex'] = impulsLoudWZ05Ex
            dataByStim.loc[renderNames[ii], 'ImpulsLoudWZPowAvg'] = impulsLoudWZPowAvg

        elif filenames[ii].find("ImpulsLoudWECMA") != -1:
            filedata = io.loadmat(file)['ImpulsLoudWECMATDep']
            
            # Calculate Willemen et al overall impulsive loudness from
            # time-dependent values
            impulsLoudWECMATDep = pd.DataFrame(filedata[1:, 1:],
                                                  index=filedata[1:, 0],
                                                  columns=["Monaural"])
            # mask for start/end skip
            impulsLoudWECMATDepMask = impulsLoudWECMATDep.loc[(impulsLoudWECMATDep.index.values
                                                               > start_skipT)
                                                              & (impulsLoudWECMATDep.index.values
                                                                 < impulsLoudWECMATDep.index.values.max()
                                                                 - end_skipT), :]

            # overall (averaged) impulsive loudness
            impulsLoudWECMAAvg = impulsLoudWECMATDepMask.mean().iloc[0]

            # overall (95th percentile = 5% exceeded) impulsive loudness
            impulsLoudWECMA05Ex = np.percentile(impulsLoudWECMATDepMask, q=95, axis=0,
                                                method='median_unbiased')

            # overall (power-averaged) impulsive loudness
            impulsLoudWECMAPowAvg = (impulsLoudWECMATDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'ImpulsLoudWECMAAvg'] = impulsLoudWECMAAvg
            dataByStim.loc[renderNames[ii], 'ImpulsLoudWECMA05Ex'] = impulsLoudWECMA05Ex
            dataByStim.loc[renderNames[ii], 'ImpulsLoudWECMAPowAvg'] = impulsLoudWECMAPowAvg

        # end of if section for absolute SQMs

        elif filenames[ii].find("PartLoudSHM") != -1:
            filedata = io.loadmat(file)['PartLoudSHMSpecTDep']

            # Calculate Sottek Hearing Model partial loudness from
            # specific time-dependent partial loudness
            specPartLoudSHM = pd.DataFrame(filedata[1:, 1:],
                                           index=filedata[1:, 0],
                                           columns=filedata[0, 1:])

            # time-dependent partial loudness (ECMA-418-2:2025 Equation 116)
            partLoudSHMTDep = specPartLoudSHM.sum(axis=1)*bandDiff0p5

            # mask for start/end skip
            partLoudSHMTDepMask = partLoudSHMTDep.loc[(partLoudSHMTDep.index.values
                                                       > start_skipT)
                                                      & (partLoudSHMTDep.index.values
                                                         < partLoudSHMTDep.index.values.max()
                                                         - end_skipT)]

            # overall (power-averaged) partial loudness (ECMA-418-2:2022 Equation 117)
            partLoudSHMPowAvg = (partLoudSHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'PartLoudSHMPowAvg'] = partLoudSHMPowAvg

        
        elif filenames[ii].find("PartTonLdSHM") != -1:
            filedata = io.loadmat(file)['PartTonLdSHMSpecTDep']

            # Calculate Sottek Hearing Model partial tonal loudness from
            # specific time-dependent partial tonal loudness
            specPartTonLdSHM = pd.DataFrame(filedata[1:, 1:],
                                             index=filedata[1:, 0],
                                             columns=filedata[0, 1:])

            # time-dependent partial tonal loudness (ECMA-418-2:2025 Equation 116)
            partTonLdSHMTDep = specPartTonLdSHM.sum(axis=1)*bandDiff0p5

            # mask for start/end skip
            partTonLdSHMTDepMask = partTonLdSHMTDep.loc[(partTonLdSHMTDep.index.values
                                                         > start_skipT)
                                                        & (partTonLdSHMTDep.index.values
                                                           < partTonLdSHMTDep.index.values.max()
                                                           - end_skipT)]

            # overall (power-averaged) partial tonal loudness (ECMA-418-2:2022 Equation 117)
            partTonLdSHMPowAvg = (partTonLdSHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'PartTonLdSHMPowAvg'] = partTonLdSHMPowAvg


        # calculation section for SQM differences
        # exclude FluctFZ and LoudECMA
        if (filenames[ii].find("FluctFZ") == -1) and (filenames[ii].find("LoudECMA") == -1): 

            # NOTE: THIS SECTION RELIES ON THE ALPHABETIC ORDER OF THE STIMULI FILES AS
            # ORIGINALLY NAMED: EACH AMBIENT SOUND FILE PRECEDING THE CORRESPONDING
            # COMBINED STIMULI FILES. IF THE ORDERING OR FILE NAMING IS CHANGED, THIS
            # CALCULATION WILL BE INVALIDATED AND *MAY ALSO* CAUSE AN ERROR.
            # a rolling 50 ms window is applied to average the SQM values over time -
            # this is to reduce uncertainty due to imperfect time-alignment between the
            # ambient vs combined stimuli files (all are from recordings, so there will
            # be some slippage due to imperfect editing)
            if renderNames[ii].__contains__("_Baseline_"):
                # NOTE: we could dropna() the first <windowT values, but these will be
                # ignored anyway in the statistical analysis, assuming start_skipT >
                # windowT

                # calculate moving average values for ambient stimulus
                if filenames[ii].find("RoughFZ") != -1:
                    ambSpecRoughFZMovAvg = specRoughFZ.rolling(window=int(np.ceil(sampleRateRoughFZ*windowT))).mean()

                elif filenames[ii].find("TonalECMA") != -1:
                    ambSpecTonalECMAMovAvg = specTonalECMA.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean()

                    # tonal annoyance weighted
                    ambSpecTonalAwSHMMovAvg = specTonalAwSHM.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean()

                elif filenames[ii].find("TonLdECMA") != -1:
                    ambSpecTonLdECMAMovAvg = specTonLdECMA.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()

                elif filenames[ii].find("RoughECMA") != -1:
                    ambSpecRoughECMAMovAvg = specRoughECMA.rolling(window=int(np.ceil(sampleRateModECMA*windowT))).mean()

                elif filenames[ii].find("FluctOV") != -1:
                    ambSpecFluctOVMovAvg = specFluctOV.rolling(window=int(np.ceil(sampleRateFluctOV*windowT))).mean()

                elif filenames[ii].find("ImpulsLoudWZ") != -1:
                    ambImpulsLoudWZTDepMovAvg = impulsLoudWZTDep.rolling(window=int(np.ceil(sampleRateLoudISO1HighRes*windowT))).mean()
                
                elif filenames[ii].find("ImpulsLoudWECMA") != -1:
                    ambImpulsLoudWECMATDepMovAvg = impulsLoudWECMATDep.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()

            elif renderNames[ii].__contains__("Street") or renderNames[ii].__contains__("Park"):
                # calculate moving average values for combined stimulus

                if filenames[ii].find("RoughFZ") != -1:
                    specRoughFZMovAvg = specRoughFZ.rolling(window=int(np.ceil(sampleRateRoughFZ*windowT))).mean()

                    # calculate differences and make negative values 0
                    dSpecRoughFZ = np.maximum(specRoughFZMovAvg
                                              - ambSpecRoughFZMovAvg, 0)

                    # calculate aggregated difference values

                    # time-dependent roughness
                    dRoughFZTDep = bandDiff0p5*dSpecRoughFZ.sum(axis=1).to_frame(name="Monaural")

                    # mask for start/end skip
                    dRoughFZTDepMask = dRoughFZTDep.loc[(dRoughFZTDep.index.values
                                                         > start_skipT)
                                                        & (dRoughFZTDep.index.values
                                                           < dRoughFZTDep.index.values.max()
                                                           - end_skipT)]

                    # overall (90th percentile = 10% exceeded) roughness
                    dRoughFZ10Ex = np.percentile(dRoughFZTDepMask, q=90, axis=0,
                                                 method='median_unbiased')

                    # overall (95th percentile = 5% exceeded) roughness
                    dRoughFZ05Ex = np.percentile(dRoughFZTDepMask, q=95, axis=0,
                                                 method='median_unbiased')

                    # add results to output DataFrame
                    dataByStim.loc[renderNames[ii], 'dRoughFZ10Ex'] = dRoughFZ10Ex
                    dataByStim.loc[renderNames[ii], 'dRoughFZ05Ex'] = dRoughFZ05Ex

                elif filenames[ii].find("TonalECMA") != -1:
                    specTonalECMAMovAvg = specTonalECMA.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean()
                    
                    # calculate differences and make negative values 0
                    dSpecTonalECMA = np.maximum(specTonalECMAMovAvg
                                                - ambSpecTonalECMAMovAvg, 0)

                    # time-dependent tonality (max, not integration)
                    dTonalECMATDep = dSpecTonalECMAL.max(axis=1).to_frame(name="Monaural")
                    
                    # time-dependent integrated tonality
                    dTonalSHMIntTDep = (dSpecTonalECMA.sum(axis=1)*bandDiff0p5*0.348088948583815).to_frame(name="Monaural")

                    # mask for start/end skip and values <= 0.02
                    dTonalECMATDepMask = dTonalECMATDep.loc[(dTonalECMATDep.index.values
                                                             > start_skipT)
                                                            & (dTonalECMATDep.index.values
                                                               < dTonalECMATDep.index.values.max()
                                                               - end_skipT)
                                                            & (dTonalECMATDep.loc[:, 0].values
                                                               > 0.02), 0]

                    # mask for start/end skip and values <= 0.02
                    # NOTE: uses tonality magnitude mask from ECMA tonality
                    dTonalSHMIntTDepMask = dTonalSHMIntTDep.loc[(dTonalSHMIntTDep.index.values
                                                                 > start_skipT)
                                                                & (dTonalSHMIntTDep.index.values
                                                                   < dTonalSHMIntTDep.index.values.max()
                                                                   - end_skipT)
                                                                & (dTonalECMATDep.loc[:, 0].values
                                                                   > 0.02), 0]  # see NOTE above

                    # time-averaged tonality (omitting T<=0.02)
                    dTonalECMAAvg = dTonalECMATDepMaskL.mean(axis=0).iloc[0]

                    # 5% exceeded tonality (automatically omitting T<=0.02)
                    dTonalECMA05Ex = np.percentile(dTonalECMATDepMaskL, q=95, axis=0,
                                                   method='median_unbiased')

                    # time-averaged integrated tonality (omitting T<=0.02)
                    dTonalSHMIntAvg = dTonalSHMIntTDepMask.mean(axis=0).iloc[0]

                    # 5% exceeded integated tonality (automatically omitting T<=0.02)
                    dTonalSHMInt05Ex = np.percentile(dTonalSHMIntTDepMask, q=95, axis=0,
                                                     method='median_unbiased')

                    # add results to output DataFrame
                    dataByStim.loc[renderNames[ii], 'dTonalECMAAvg'] = dTonalECMAAvg
                    dataByStim.loc[renderNames[ii], 'dTonalECMA05Ex'] = dTonalECMA05Ex
                    dataByStim.loc[renderNames[ii], 'dTonalSHMIntAvg'] = dTonalSHMIntAvg
                    dataByStim.loc[renderNames[ii], 'dTonalSHMInt05Ex'] = dTonalSHMInt05Ex

                    # tonal annoyance-weighted
                    specTonalAwSHMMovAvg = specTonalAwSHM.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean()
        
                    # calculate differences and make negative values 0
                    dSpecTonalAwSHM = np.maximum(specTonalAwSHMMovAvg
                                                 - ambSpecTonalAwSHMLMovAvg, 0)

                    # time-dependent tonal annoyance-weighted tonality
                    dTonalAwSHMTDep = dSpecTonalAwSHML.max(axis=1).to_frame(name="Monaural")

                    # time-dependent integrated tonal annoyance-weighted tonality
                    dTonalAwSHMIntTDep = (dSpecTonalAwSHM.sum(axis=1)*bandDiff0p5*0.348088948583815).to_frame(name="Monaural")
                    
                    # mask for start/end skip and values <= 0.02
                    dTonalAwSHMTDepMask = dTonalAwSHMTDep.loc[(dTonalAwSHMTDep.index.values
                                                               > start_skipT)
                                                              & (dTonalAwSHMTDep.index.values
                                                                 < dTonalAwSHMTDep.index.values.max()
                                                                 - end_skipT)
                                                              & (dTonalAwSHMTDep.loc[:, 0].values
                                                                 > 0.02), 0]

                    # mask for start/end skip and values <= 0.02
                    # NOTE: uses tonality magnitude mask from ECMA tonality
                    dTonalAwSHMIntTDepMask = dTonalAwSHMIntTDep.loc[(dTonalAwSHMIntTDep.index.values
                                                                     > start_skipT)
                                                                    & (dTonalAwSHMIntTDep.index.values
                                                                       < dTonalAwSHMIntTDep.index.values.max()
                                                                       - end_skipT)
                                                                    & (dTonalAwSHMTDep.loc[:, 0].values
                                                                       > 0.02), 0]  # see NOTE above

                    # time-averaged tonality (omitting T<=0.02)
                    dTonalAwSHMAvg = dTonalAwSHMTDepMask.mean(axis=0).iloc[0]

                    # 5% exceeded tonality (automatically omitting T<=0.02)
                    dTonalAwSHM05ExMaxLR = np.percentile(dTonalAwSHMTDepMask, q=95, axis=0,
                                                         method='median_unbiased')

                    # time-averaged integrated tonality (omitting T<=0.02)
                    dTonalAwSHMIntAvgMaxLR = dTonalAwSHMIntTDepMaskL.mean(axis=0).iloc[0]

                    # 5% exceeded integrated tonality (automatically omitting T<=0.02)
                    dTonalAwSHMInt05Ex = np.percentile(dTonalAwSHMIntTDepMaskL, q=95, axis=0,
                                                       method='median_unbiased')

                    # add results to output DataFrame
                    dataByStim.loc[renderNames[ii], 'dTonalAwSHMAvg'] = dTonalAwSHMAvg
                    dataByStim.loc[renderNames[ii], 'dTonalAwSHM05Ex'] = dTonalAwSHM05Ex
                    dataByStim.loc[renderNames[ii], 'dTonalAwSHMIntAvg'] = dTonalAwSHMIntAvg
                    dataByStim.loc[renderNames[ii], 'dTonalAwSHMInt05Ex'] = dTonalAwSHMInt05Ex
                    
                elif filenames[ii].find("TonLdECMA") != -1:
                    specTonLdECMAMovAvg = specTonLdECMA.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()
                    
                    # calculate differences and make negative values 0
                    dSpecTonLdECMA = np.maximum(specTonLdECMAMovAvg
                                                - ambSpecTonLdECMAMovAvg, 0)

                    # time-dependent tonal loudness (ECMA-418-2:2025 Equation 116)
                    dTonLdECMATDep = dSpecTonLdECMA.sum(axis=1)*bandDiff0p5

                    # mask for start/end skip
                    dTonLdECMATDepMask = dTonLdECMATDep.loc[(dTonLdECMATDep.index.values
                                                             > start_skipT)
                                                            & (dTonLdECMATDep.index.values
                                                               < dTonLdECMATDep.index.values.max()
                                                               - end_skipT)]

                    # overall (power-averaged) tonal loudness (ECMA-418-2:2025
                    # Equation 117)
                    dTonLdECMAPowAvg = dTonLdECMATDepMask.pow(1/np.log10(2)).mean()**np.log10(2)

                    # 5% exceeded tonal loudness
                    dTonLdECMA05Ex = np.percentile(dTonLdECMATDepMask, q=95, axis=0,
                                                   method='median_unbiased')

                    # add results to output DataFrame
                    dataByStim.loc[renderNames[ii], 'dTonLdECMAPowAvg'] = dTonLdECMAPowAvg
                    dataByStim.loc[renderNames[ii], 'dTonLdECMA05Ex'] = dTonLdECMA05Ex

                elif filenames[ii].find("RoughECMA") != -1:
                    specRoughECMAMovAvg = specRoughECMA.rolling(window=int(np.ceil(sampleRateModECMA*windowT))).mean()
        
                    # calculate differences and make negative values 0
                    dSpecRoughECMA = np.maximum(specRoughECMAMovAvg
                                                - ambSpecRoughECMAMovAvg, 0)

                    # binaural time-dependent roughness
                    dRoughECMATDep = dSpecRoughECMA.sum(axis=1)*bandDiff0p5

                    # mask for start/end skip
                    dRoughECMATDepMask = dRoughECMATDep.loc[(dRoughECMATDep.index.values
                                                             > start_skipT)
                                                            & (dRoughECMATDep.index.values
                                                               < dRoughECMATDep.index.values.max()
                                                               - end_skipT)]

                    # overall (90th percentile = 10% exceeded) roughness
                    dRoughECMA10Ex = np.percentile(dRoughECMATDepMask, q=90, axis=0,
                                                   method='median_unbiased')
                    
                    # overall (95th percentile = 5% exceeded) roughness
                    dRoughECMA05Ex = np.percentile(dRoughECMATDepMask, q=95, axis=0,
                                                   method='median_unbiased')

                    # add results to output DataFrame
                    dataByStim.loc[renderNames[ii], 'dRoughECMA10Ex'] = dRoughECMA10Ex
                    dataByStim.loc[renderNames[ii], 'dRoughECMA05Ex'] = dRoughECMA05Ex

                elif filenames[ii].find("FluctOV") != -1:
                    specFluctOVMovAvg = specFluctOV.rolling(window=int(np.ceil(sampleRateFluctOV*windowT))).mean()

                    # calculate differences and make negative values 0
                    dSpecFluctOV = np.maximum(specFluctOVMovAvg
                                              - ambSpecFluctOVMovAvg, 0)

                    # time-dependent integrated fluctuation strength
                    dFluctOVTDep = dSpecFluctOVL.sum(axis=1)*bandDiff0p5
                    
                    # fluctuation strength masked for start/end skip
                    dFluctOVTDepMask = dFluctOVTDep.loc[(dFluctOVTDep.index.values
                                                         > start_skipT).transpose()
                                                        & (dFluctOVTDep.index.values
                                                           < dFluctOVTDep.index.values.max()
                                                           - end_skipT).transpose()]

                    # overall 5% exceeded fluctuation strength 
                    dFluctOV05Ex = np.percentile(dFluctOVTDepMask, q=95, axis=0,
                                                 method='median_unbiased')

                    # overall 10% exceeded fluctuation strength 
                    dFluctOV10Ex = np.percentile(dFluctOVTDepMask, q=90, axis=0,
                                                 method='median_unbiased')
                    
                    # add results to output DataFrame
                    dataByStim.loc[renderNames[ii], 'dFluctOV10Ex'] = dFluctOV10Ex
                    dataByStim.loc[renderNames[ii], 'dFluctOV05Ex'] = dFluctOV05Ex
                
                elif filenames[ii].find("ImpulsLoudWZ") != -1:
                    impulsLoudWZTDepMovAvg = impulsLoudWZTDep.rolling(window=int(np.ceil(sampleRateLoudISO1HighRes*windowT))).mean()
                    
                    # calculate differences and make negative values 0
                    dImpulsLoudWZTDep = np.maximum(impulsLoudWZTDepMovAvg
                                                   - ambImpulsLoudWZTDepMovAvg, 0)
                    
                    # mask for start/end skip
                    dImpulsLoudWZTDepMask = dImpulsLoudWZTDep.loc[(dImpulsLoudWZTDep.index.values
                                                                   > start_skipT)
                                                                  & (dImpulsLoudWZTDep.index.values
                                                                     < dImpulsLoudWZTDep.index.values.max()
                                                                     - end_skipT), :]

                    # overall (averaged) impulsive loudness
                    dImpulsLoudWZAvg = dImpulsLoudWZTDepMask.mean(axis=0).iloc[0]
                    
                    # overall (95th percentile = 5% exceeded) impulsive loudness
                    dImpulsLoudWZ05Ex = np.percentile(dImpulsLoudWZTDepMask, q=95, axis=0,
                                                      method='median_unbiased')
                    
                    # overall (power-averaged) impulsive loudness
                    dImpulsLoudWZPowAvg = (dImpulsLoudWZTDepMask.pow(1/np.log10(2)).mean(axis=0)**np.log10(2)).iloc[0]

                    # add results to output DataFrame
                    dataByStim.loc[renderNames[ii], 'dImpulsLoudWZAvg'] = dImpulsLoudWZAvg
                    dataByStim.loc[renderNames[ii], 'dImpulsLoudWZ05Ex'] = dImpulsLoudWZ05Ex
                    dataByStim.loc[renderNames[ii], 'dImpulsLoudWZPowAvg'] = dImpulsLoudWZPowAvg

                elif filenames[ii].find("ImpulsLoudWECMA") != -1:
                    impulsLoudWECMATDepMovAvg = impulsLoudWECMATDep.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()
                    
                    # calculate differences and make negative values 0
                    dImpulsLoudWECMATDep = np.maximum(impulsLoudWECMATDepMovAvg
                                                      - ambImpulsLoudWECMATDepMovAvg, 0)
                    
                    # mask for start/end skip
                    dImpulsLoudWECMATDepMask = dImpulsLoudWECMATDep.loc[(dImpulsLoudWECMATDep.index.values
                                                                         > start_skipT)
                                                                        & (dImpulsLoudWECMATDep.index.values
                                                                           < dImpulsLoudWECMATDep.index.values.max()
                                                                           - end_skipT), :]

                    # overall (averaged) impulsive loudness
                    dImpulsLoudWECMAAvg = dImpulsLoudWECMATDepMask.mean().iloc[0]
                    
                    # overall (95th percentile = 5% exceeded) impulsive loudness
                    dImpulsLoudWECMA05Ex = np.percentile(dImpulsLoudWECMATDepMask, q=95, axis=0,
                                                         method='median_unbiased')
                    
                    # overall (power-averaged) impulsive loudness
                    dImpulsLoudWECMAPowAvg = (dImpulsLoudWECMATDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

                    # add results to output DataFrame
                    dataByStim.loc[renderNames[ii], 'dImpulsLoudWECMAAvg'] = dImpulsLoudWECMAAvg
                    dataByStim.loc[renderNames[ii], 'dImpulsLoudWECMA05Ex'] = dImpulsLoudWECMA05Ex
                    dataByStim.loc[renderNames[ii], 'dImpulsLoudWECMAPowAvg'] = dImpulsLoudWECMAPowAvg

        # end of if section for SQM differences
    # end of if branch for .mat files

    elif file[-4:] == "xlsx":
        workbookdata = pd.read_excel(io=file, sheet_name=None, engine='calamine')

        # Calculate quasi-Zwicker (ISO 532-1) overall loudness from 
        # time-dependent loudness
        loudQZ5321TDep = pd.DataFrame(workbookdata['LoudQZ5321'].iloc[0:, 1:].values,
                                      columns=workbookdata['LoudQZ5321'].iloc[0, 1:].index,
                                      index=workbookdata['LoudQZ5321'].iloc[:, 0])

        # mask for start/end skip
        loudQZ5321TDepMask = loudQZ5321TDep.loc[(loudQZ5321TDep.index.values
                                               > start_skipT)
                                              & (loudQZ5321TDep.index.values
                                                 < loudQZ5321TDep.index.values.max()
                                                 - end_skipT), :]

        # overall (power-averaged) loudness
        loudQZ5321PowAvg = (loudQZ5321TDepMask.pow(1/np.log10(2)).mean(axis=0)**np.log10(2)).iloc[0]

        # overall 5% exceeded loudness
        loudQZ532105Ex = np.percentile(loudQZ5321TDepMask, q=95, axis=0,
                                       method='median_unbiased')
        
        # Calculate overall Aures+quasi-Zwicker (ISO 532-3 adjusted) loudness from
        # time-dependent loudness
        loudQZ5323TDep = pd.DataFrame(workbookdata['LoudQZ5323'].iloc[0:, 1:].values,
                                     columns=workbookdata['LoudQZ5323'].iloc[0, 1:].index,
                                     index=workbookdata['LoudQZ5323'].iloc[:, 0])

        # mask for start/end skip
        loudQZ5323TDepMask = loudQZ5323TDep.loc[(loudQZ5323TDep.index.values
                                                 > start_skipT)
                                                & (loudQZ5323TDep.index.values
                                                   < loudQZ5323TDep.index.values.max()
                                                   - end_skipT), :]

        # overall (power-averaged) loudness
        loudQZ5323PowAvg = (loudQZ5323TDepMask.pow(1/np.log10(2)).mean(axis=0)**np.log10(2)).iloc[0]

        # overall 5% exceeded loudness
        loudQZ532305Ex = np.percentile(loudQZ5323TDepMask, q=95, axis=0,
                                       method='median_unbiased')

        # Calculate overall Aures+quasi-Zwicker (ECMA-418-2 ear filter and
        # transformation) loudness from time-dependent loudness
        loudQZ4182TDep = pd.DataFrame(workbookdata['LoudQZ4182'].iloc[0:, 1:].values,
                                      columns=workbookdata['LoudQZ4182'].iloc[0, 1:].index,
                                      index=workbookdata['LoudQZ4182'].iloc[:, 0])

        # mask for start/end skip
        loudQZ4182TDepMask = loudQZ4182TDep.loc[(loudQZ4182TDep.index.values
                                                 > start_skipT)
                                                & (loudQZ4182TDep.index.values
                                                   < loudQZ4182TDep.index.values.max()
                                                   - end_skipT), :]
        
        # overall (power-averaged) loudness
        loudQZ4182PowAvg = (loudQZ4182TDepMask.pow(1/np.log10(2)).mean(axis=0)**np.log10(2)).iloc[0]
    
        # overall 5% exceeded loudness
        loudQZ418205Ex = np.percentile(loudQZ4182TDepMask, q=95, axis=0,
                                       method='median_unbiased')

        # Calculate overall Aures+quasi-Zwicker sharpness from
        # time-dependent sharpness
        sharpAQZ5321TDep = pd.DataFrame(workbookdata['SharpAuresQZ5321'].iloc[0:, 1:].values,
                                        columns=workbookdata['SharpAuresQZ5321'].iloc[0, 1:].index,
                                        index=workbookdata['SharpAuresQZ5321'].iloc[:, 0])

        # mask for start/end skip
        sharpAQZ5321TDepMask = sharpAQZ5321TDep.loc[(sharpAQZ5321TDep.index.values
                                                     > start_skipT)
                                                    & (sharpAQZ5321TDep.index.values
                                                       < sharpAQZ5321TDep.index.values.max()
                                                       - end_skipT)]

        # overall (power-averaged) sharpness
        sharpAQZ5321PowAvg = (sharpAQZ5321TDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

        # overall 5% exceeded sharpness
        sharpAQZ532105Ex = np.percentile(sharpAQZ5321TDepMask, q=95, axis=0,
                                         method='median_unbiased')

        # Calculate overall Aures+quasi-Zwicker (ISO 532-3 adjusted) sharpness from
        # time-dependent sharpness
        sharpAQZ5323TDep = pd.DataFrame(workbookdata['SharpAuresQZ5323'].iloc[0:, 1:].values,
                                        columns=workbookdata['SharpAuresQZ5323'].iloc[0, 1:].index,
                                        index=workbookdata['SharpAuresQZ5323'].iloc[:, 0])

        # mask for start/end skip
        sharpAQZ5323TDepMask = sharpAQZ5323TDep.loc[(sharpAQZ5323TDep.index.values
                                                     > start_skipT)
                                                    & (sharpAQZ5323TDep.index.values
                                                       < sharpAQZ5323TDep.index.values.max()
                                                       - end_skipT)]

        # overall (power-averaged) sharpness
        sharpAQZ5323PowAvg = (sharpAQZ5323TDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

        # overall 5% exceeded sharpness
        sharpAQZ532305Ex = np.percentile(sharpAQZ5323TDepMask, q=95, axis=0,
                                         method='median_unbiased')
        
        # Calculate overall Aures+quasi-Zwicker (ECMA-418-2 ear filter and
        # transformation) sharpness from time-dependent sharpness
        sharpAQZ4182TDep = pd.DataFrame(workbookdata['SharpAuresQZ4182'].iloc[0:, 1:].values,
                                        columns=workbookdata['SharpAuresQZ4182'].iloc[0, 1:].index,
                                        index=workbookdata['SharpAuresQZ4182'].iloc[:, 0])

        # mask for start/end skip
        sharpAQZ4182TDepMask = sharpAQZ4182TDep.loc[(sharpAQZ4182TDep.index.values
                                                     > start_skipT)
                                                    & (sharpAQZ4182TDep.index.values
                                                       < sharpAQZ4182TDep.index.values.max()
                                                       - end_skipT)]

        # overall (power-averaged) sharpness
        sharpAQZ4182PowAvg = (sharpAQZ4182TDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

        # overall 5% exceeded sharpness
        sharpAQZ418205Ex = np.percentile(sharpAQZ4182TDepMask, q=95, axis=0,
                                         method='median_unbiased')
        
        # Calculate overall Aures tonality from 2-channel time-dependent
        # tonality
        tonalAurTDep = pd.DataFrame(workbookdata['TonalAures'].iloc[0:, 1:].values,
                                    columns=workbookdata['TonalAures'].iloc[0, 1:].index,
                                    index=workbookdata['TonalAures'].iloc[:, 0])

        # mask for start/end skip
        tonalAurTDepMask = tonalAurTDep.loc[(tonalAurTDep.index.values
                                             > start_skipT)
                                            & (tonalAurTDep.index.values
                                               < tonalAurTDep.index.values.max()
                                               - end_skipT), :]

        # overall mean tonality
        tonalAurAvg = tonalAurTDepMask.mean(axis=0).iloc[0]

        # overall 5% exceeded tonality
        tonalAur05Ex = np.percentile(tonalAurTDepMask, q=95, axis=0,
                                     method='median_unbiased')

        # overall 10% exceeded tonality
        tonalAur10Ex = np.percentile(tonalAurTDepMask, q=90, axis=0,
                                     method='median_unbiased')

        # Calculate Daniel & Weber overall roughness from specific roughness
        specRoughDW = pd.DataFrame(workbookdata['RoughDanWeb'].iloc[:, 1:].values.T,
                                   columns=workbookdata['RoughDanWeb'].iloc[:, 0],
                                   index=workbookdata['RoughDanWeb'].iloc[0, 1:].index)

        # time-dependent roughness
        roughDWTDep = bandDiff0p5*specRoughDW.sum(axis=0)

        # mask for start/end skip
        roughDWTDepMask = roughDWTDep.loc[(roughDWTDep.index.values
                                            > start_skipT)
                                          & (roughDWTDep.index.values
                                              < roughDWTDep.index.values.max()
                                              - end_skipT)]

        # overall (90th percentile = 10% exceeded) roughness
        roughDW10ExMaxLR = np.percentile(roughDWTDepMask, q=90, axis=0,
                                         method='median_unbiased')

        # overall (95th percentile = 5% exceeded) roughness
        roughDW05ExMaxLR = np.percentile(roughDWTDepMask, q=95, axis=0,
                                         method='median_unbiased')

        # Calculate overall Aures+Sottek Hearing Model sharpness from
        # time-dependent sharpness
        sharpASHMTDep = pd.DataFrame(workbookdata['SharpAuresSHM'].iloc[0:, 1:].values,
                                     columns=workbookdata['SharpAuresSHM'].iloc[0, 1:].index,
                                     index=workbookdata['SharpAuresSHM'].iloc[:, 0])

        # mask for start/end skip
        sharpASHMTDepMask = sharpASHMTDep.loc[(sharpASHMTDep.index.values
                                               > start_skipT)
                                              & (sharpASHMTDep.index.values
                                                 < sharpASHMTDep.index.values.max()
                                                 - end_skipT)]

        # overall (power-averaged) sharpness
        sharpASHMPowAvg = sharpASHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)
        sharpASHMPowAvg = sharpASHMPowAvg.iloc[0]

        # overall 5% exceeded sharpness
        sharpASHM05Ex = np.percentile(sharpASHMTDepMask, q=95, axis=0,
                                      method='median_unbiased')

        # Calculate overall Aures+Sottek Hearing Model tonal sharpness from
        # time-dependent tonal sharpness
        tonShpASHMTDep = pd.DataFrame(workbookdata['TonalSharpAuresSHM'].iloc[0:, 1:].values,
                                      columns=workbookdata['TonalSharpAuresSHM'].iloc[0, 1:].index,
                                      index=workbookdata['TonalSharpAuresSHM'].iloc[:, 0])

        # mask for start/end skip
        tonShpASHMTDepMask = tonShpASHMTDep.loc[(tonShpASHMTDep.index.values
                                                 > start_skipT)
                                                & (tonShpASHMTDep.index.values
                                                   < tonShpASHMTDep.index.values.max()
                                                   - end_skipT)]

        # overall (power-averaged) tonal sharpness
        tonShpASHMPowAvg = (tonShpASHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

        # overall 5% exceeded tonal sharpness
        tonShpASHM05Ex = np.percentile(tonShpASHMTDepMask, q=95, axis=0,
                                       method='median_unbiased')

        # add results to output DataFrame
        dataByStim.loc[renderNames[ii], 'LoudQZ5321PowAvg'] = loudQZ5321PowAvg
        dataByStim.loc[renderNames[ii], 'LoudQZ532105Ex'] = loudQZ532105Ex
        dataByStim.loc[renderNames[ii], 'LoudQZ5323PowAvg'] = loudQZ5323PowAvg
        dataByStim.loc[renderNames[ii], 'LoudQZ532305Ex'] = loudQZ532305Ex
        dataByStim.loc[renderNames[ii], 'LoudQZ5323PowAvg'] = loudQZ5323PowAvg
        dataByStim.loc[renderNames[ii], 'LoudQZ4182PowAvg'] = loudQZ4182PowAvg
        dataByStim.loc[renderNames[ii], 'LoudQZ418205Ex'] = loudQZ418205Ex
        dataByStim.loc[renderNames[ii], 'LoudQZ4182PowAvg'] = loudQZ4182PowAvg
        dataByStim.loc[renderNames[ii], 'SharpAurQZ5321PowAvg'] = sharpAQZ5321PowAvg
        dataByStim.loc[renderNames[ii], 'SharpAurQZ532105Ex'] = sharpAQZ532105Ex
        dataByStim.loc[renderNames[ii], 'SharpAurQZ5323PowAvg'] = sharpAQZ5323PowAvg
        dataByStim.loc[renderNames[ii], 'SharpAurQZ532305Ex'] = sharpAQZ532305Ex
        dataByStim.loc[renderNames[ii], 'SharpAurQZ4182PowAvg'] = sharpAQZ4182PowAvg
        dataByStim.loc[renderNames[ii], 'SharpAurQZ418205Ex'] = sharpAQZ418205Ex
        dataByStim.loc[renderNames[ii], 'TonalAur05Ex'] = tonalAur05Ex
        dataByStim.loc[renderNames[ii], 'TonalAur10Ex'] = tonalAur10Ex
        dataByStim.loc[renderNames[ii], 'TonalAurAvg'] = tonalAurAvg
        dataByStim.loc[renderNames[ii], 'RoughDW10Ex'] = roughDW10Ex
        dataByStim.loc[renderNames[ii], 'RoughDW05Ex'] = roughDW05Ex
        dataByStim.loc[renderNames[ii], 'SharpAurSHMPowAvg'] = sharpASHMPowAvg
        dataByStim.loc[renderNames[ii], 'SharpAurSHM05Ex'] = sharpASHM05Ex
        dataByStim.loc[renderNames[ii], 'TonShpAurSHMPowAvg'] = tonShpASHMPowAvg
        dataByStim.loc[renderNames[ii], 'TonShpAurSHM05Ex'] = tonShpASHM05Ex

        # calculation section for Aures+Sottek Hearing Model partial tonal sharpness        
        if workbookdata.keys().__contains__('PartTonShpAuresSHM'):

            # Calculate overall Aures+Sottek Hearing Model partial tonal sharpness from
            # time-dependent partial tonal sharpness
            partTonShpASHMTDep = pd.DataFrame(workbookdata['PartTonShpAuresSHM'].iloc[0:, 1:].values,
                                              columns=workbookdata['PartTonShpAuresSHM'].iloc[0, 1:].index,
                                              index=workbookdata['PartTonShpAuresSHM'].iloc[:, 0])

            # mask for start/end skip
            partTonShpASHMTDepMask = partTonShpASHMTDep.loc[(partTonShpASHMTDep.index.values
                                                             > start_skipT)
                                                            & (partTonShpASHMTDep.index.values
                                                               < partTonShpASHMTDep.index.values.max()
                                                               - end_skipT)]

            # overall (power-averaged) partial tonal sharpness
            partTonShpASHMPowAvg = partTonShpASHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)
            partTonShpASHMPowAvgBin = partTonShpASHMPowAvgBin.iloc[0]

            # overall 5% exceeded partial tonal sharpness
            partTonShpASHM05Ex = np.percentile(partTonShpASHMTDepMask, q=95, axis=0,
                                               method='median_unbiased')

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'PartTonShpAurSHMPowAvg'] = partTonShpASHMPowAvg
            dataByStim.loc[renderNames[ii], 'PartTonShpAurSHM05Ex'] = partTonShpASHM05Ex

        # calculation section for SQM differences
        # NOTE: THIS SECTION RELIES ON THE ALPHABETIC ORDER OF THE STIMULI FILES AS
        # ORIGINALLY NAMED: EACH AMBIENT SOUND FILE PRECEDING THE CORRESPONDING
        # COMBINED STIMULI FILES. IF THE ORDERING OR FILE NAMING IS CHANGED, THIS
        # CALCULATION WILL BE INVALIDATED AND *MAY ALSO* CAUSE AN ERROR.
        # a rolling 50 ms window is applied to average the SQM values over time -
        # this is to reduce uncertainty due to imperfect time-alignment between the
        # ambient vs combined stimuli files (all are from recordings, so there will
        # be some slippage due to imperfect editing)
        if renderNames[ii].__contains__("_Baseline_"):
            # NOTE: we could dropna() the first <windowT values, but these will be
            # ignored anyway in the statistical analysis, assuming start_skipT >
            # windowT

            # calculate moving average values for ambient stimulus
            ambSharpASHMTDepMovAvg = sharpASHMTDep.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()
            ambTonShpASHMTDepMovAvg = tonShpASHMTDep.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()

        elif renderNames[ii].__contains__("Street") or renderNames[ii].__contains__("Park"):
            # calculate moving average values for combined stimulus
            sharpASHMTDepMovAvg = sharpASHMTDep.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()
            tonShpASHMTDepMovAvg = tonShpASHMTDep.rolling(window=int(np.ceil(sampleRateLoudECMA*windowT))).mean()

            # # calculate differences and make negative values 0
            dSharpASHMTDep = np.maximum(sharpASHMTDepMovAvg
                                        - ambSharpASHMTDepMovAvg, 0)
            dTonShpASHMTDep = np.maximum(tonShpASHMTDepMovAvg
                                         - ambTonShpASHMTDepMovAvg, 0)
            # calculate aggregated difference values
            # sharpness masked for start/end skip
            dSharpASHMTDepMask = dSharpASHMTDep.loc[(dSharpASHMTDep.index.values
                                                     > start_skipT).transpose()
                                                    & (dSharpASHMTDep.index.values
                                                       < dSharpASHMTDep.index.values.max()
                                                       - end_skipT).transpose()]

            # overall (power-averaged) sharpness
            dSharpASHMPowAvg = (dSharpASHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall 5% exceeded sharpness
            dSharpASHM05Ex = np.percentile(dSharpASHMTDepMask, q=95, axis=0,
                                           method='median_unbiased')

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'dSharpAurSHMPowAvg'] = dSharpASHMPowAvg
            dataByStim.loc[renderNames[ii], 'dSharpAurSHM05Ex'] = dSharpASHM05Ex

            # tonal sharpness masked for start/end skip
            dTonShpASHMTDepMask = dTonShpASHMTDep.loc[(dTonShpASHMTDep.index.values
                                                       > start_skipT).transpose()
                                                      & (dTonShpASHMTDep.index.values
                                                         < dTonShpASHMTDep.index.values.max()
                                                         - end_skipT).transpose()]

            # overall (power-averaged) tonal sharpness
            dTonShpASHMPowAvg = (dTonShpASHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall 5% exceeded tonal sharpness
            dTonShpASHM05Ex = np.percentile(dTonShpASHMTDepMask, q=95, axis=0,
                                            method='median_unbiased')
            
            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'dTonShpAurSHMPowAvg'] = dTonShpASHMPowAvg
            dataByStim.loc[renderNames[ii], 'dTonShpAurSHM05Ex'] = dTonShpASHM05Ex

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

            if file.find("Impulsiveness") != -1:
                # Sottek Hearing Model impulsiveness
                # ----------------------------------
                filedata = np.array(pd.read_csv(file, sep="	", header=None))
                chan2row = np.arange(0, filedata.shape[0])[np.isnan(filedata[:, 0])][-1]
        
                # specific time-dependent impulsiveness
                specImpulsSHML = pd.DataFrame(filedata[1:chan2row, 1:],
                                            index=filedata[1:chan2row, 0],
                                            columns=filedata[0, 1:])
                specImpulsSHMR = pd.DataFrame(filedata[chan2row + 1:, 1:],
                                            index=filedata[chan2row + 1:, 0],
                                            columns=filedata[chan2row, 1:])
            
            elif file.find("Fluctuation") != -1:
                # ECMA-418-2:2025 Fluctuation strength
                # -----------------------------------
                filedata = np.array(pd.read_csv(file, sep="	", header=None))
                chan2row = np.arange(0, filedata.shape[0])[np.isnan(filedata[:, 0])][-1]

                # specific ECMA-418-2:2025 Fluctuation strength
                specFluctECMAL = pd.DataFrame(filedata[1:chan2row, 1:],
                                              index=filedata[1:chan2row, 0],
                                              columns=filedata[0, 1:])
                specFluctECMAR = pd.DataFrame(filedata[chan2row + 1:, 1:],
                                              index=filedata[chan2row + 1:, 0],
                                              columns=filedata[chan2row, 1:])

                # calculate binaural fluctuation strength
                specFluctECMABin = ((specFluctECMAL**2
                                     + specFluctECMAR**2)/2).pow(0.5)

                # binaural time-dependent binaural fluctuation strength
                fluctECMASHMTDepBin = specFluctECMABin.sum(axis=1)*bandDiff0p5

                # binaural fluctuation strength masked for start/end skip
                fluctECMASHMTDepBinMask = fluctECMASHMTDepBin.loc[(fluctECMASHMTDepBin.index.values
                                                                   > start_skipT).transpose()
                                                                  & (fluctECMASHMTDepBin.index.values
                                                                     < fluctECMASHMTDepBin.index.values.max()
                                                                     - end_skipT).transpose()]

                # binaural overall 10% exceeded binaural fluctuation strength
                fluctECMA10ExBin = fluctECMASHMTDepBinMask.quantile(q=0.90)

                # binaural overall 5% exceeded binaural fluctuation strength
                fluctECMA05ExBin = fluctECMASHMTDepBinMask.quantile(q=0.95)

                # add results to output DataFrame
                dataByStim.loc[renderNames[ii], 'FluctECMA10ExBin'] = fluctECMA10ExBin
                dataByStim.loc[renderNames[ii], 'FluctECMA05ExBin'] = fluctECMA05ExBin

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
                if file.find("Impulsiveness") != -1:
                    ambSpecImpulsSHMLMovAvg = specImpulsSHML.rolling(window=int(np.ceil(sampleRateImpulsSHM*windowT))).mean()
                    ambSpecImpulsSHMRMovAvg = specImpulsSHMR.rolling(window=int(np.ceil(sampleRateImpulsSHM*windowT))).mean()
                
                elif file.find("Fluctuation") != -1:
                    ambSpecFluctECMASHMLMovAvg = specFluctECMAL.rolling(window=int(np.ceil((sampleRateModECMA)*windowT))).mean()
                    ambSpecFluctECMASHMRMovAvg = specFluctECMAR.rolling(window=int(np.ceil((sampleRateModECMA)*windowT))).mean()
                
            elif renderNames[ii][0:3] in ["A1_", "A2_", "B2_"]:
                # calculate moving average values for combined stimulus
                # impulsiveness
                if file.find("Impulsiveness") != -1:
                    specImpulsSHMLMovAvg = specImpulsSHML.rolling(window=int(np.ceil(sampleRateImpulsSHM*windowT))).mean()
                    specImpulsSHMRMovAvg = specImpulsSHMR.rolling(window=int(np.ceil(sampleRateImpulsSHM*windowT))).mean()
    
                    dImpulsSHMTDep = pd.concat([np.maximum(specImpulsSHMLMovAvg
                                                - ambSpecImpulsSHMLMovAvg, 0).sum(axis=1).to_frame(name="Channel 1"),
                                                np.maximum(specImpulsSHMRMovAvg
                                                - ambSpecImpulsSHMRMovAvg, 0).sum(axis=1).to_frame(name="Channel 2")],
                                               axis=1)
        
                    # calculate aggregated difference values
                    # 2-channel impulsiveness masked for start/end skip
                    dImpulsSHMTDepMask = dImpulsSHMTDep.loc[(dImpulsSHMTDep.index.values
                                                             > start_skipT).transpose()
                                                            & (dImpulsSHMTDep.index.values
                                                               < dImpulsSHMTDep.index.values.max()
                                                               - end_skipT).transpose()]
                    
                    # max l/r overall (power-averaged) impulsiveness
                    dImpulsSHMPowAvgMaxLR = (dImpulsSHMTDepMask.pow(1/np.log10(2)).mean(axis=0)**np.log10(2)).max()

                    # max l/r overall (mean) impulsiveness
                    dImpulsSHMMeanMaxLR = dImpulsSHMTDepMask.mean(axis=0).max()
        
                    # max l/r overall 5% exceeded impulsiveness
                    dImpulsSHM05ExMaxLR = dImpulsSHMTDepMask.quantile(q=0.95, axis=0).max()
    
                    # add results to output DataFrame
                    dataByStim.loc[renderNames[ii], 'dImpulsSHMPowAvgMaxLR'] = dImpulsSHMPowAvgMaxLR
                    dataByStim.loc[renderNames[ii], 'dImpulsSHMAvgMaxLR'] = dImpulsSHMMeanMaxLR
                    dataByStim.loc[renderNames[ii], 'dImpulsSHM05ExMaxLR'] = dImpulsSHM05ExMaxLR

                # fluctuation strength
                if file.find("Fluctuation") != -1:
                    specFluctECMALMovAvg = specFluctECMAL.rolling(window=int(np.ceil((sampleRateModECMA)*windowT))).mean()
                    specFluctECMARMovAvg = specFluctECMAR.rolling(window=int(np.ceil((sampleRateModECMA)*windowT))).mean()
                
                    # calculate differences and make negative values 0
                    # fluctuation strength
                    dSpecFluctECMAL = np.maximum(specFluctECMALMovAvg
                                                 - ambSpecFluctECMASHMLMovAvg, 0)
                    dSpecFluctECMAR = np.maximum(specFluctECMARMovAvg
                                                 - ambSpecFluctECMASHMRMovAvg, 0)

                    # binaural specific fluctuation strength
                    dSpecFluctECMABin = ((dSpecFluctECMAL**2
                                         + dSpecFluctECMAR**2)/2).pow(0.5) 
                    
                    # binaural time-dependent binaural fluctuation strength
                    dFluctECMASHMTDepBin = dSpecFluctECMABin.sum(axis=1)*bandDiff0p5

                    # binaural fluctuation strength masked for start/end skip
                    dFluctECMASHMTDepBinMask = dFluctECMASHMTDepBin.loc[(dFluctECMASHMTDepBin.index.values
                                                                         > start_skipT).transpose()
                                                                        & (dFluctECMASHMTDepBin.index.values
                                                                           < dFluctECMASHMTDepBin.index.values.max()
                                                                           - end_skipT).transpose()]

                    # binaural overall 10% exceeded fluctuation strength
                    dFluctECMA10ExBin = dFluctECMASHMTDepBinMask.quantile(q=0.90)

                    # binaural overall 5% exceeded fluctuation strength
                    dFluctECMA05ExBin = dFluctECMASHMTDepBinMask.quantile(q=0.95)

                    # add results to output DataFrame
                    dataByStim.loc[renderNames[ii], 'dFluctECMA10ExBin'] = dFluctECMA10ExBin
                    dataByStim.loc[renderNames[ii], 'dFluctECMA05ExBin'] = dFluctECMA05ExBin
                    
                
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

            # max l/r overall (power-averaged) impulsiveness
            impulsSHMPowAvgMaxLR = (impulsSHMTDepMask.pow(1/np.log10(2)).mean(axis=0)**np.log10(2)).max()
 
            # max l/r overall (mean) impulsiveness
            impulsSHMMeanMaxLR = impulsSHMTDepMask.mean(axis=0).max()

            # max l/r overall 5% exceeded impulsiveness
            impulsSHM05ExMaxLR = impulsSHMTDepMask.quantile(q=0.95).max()

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'ImpulsSHMPowAvgMaxLR'] = impulsSHMPowAvgMaxLR
            dataByStim.loc[renderNames[ii], 'ImpulsSHMAvgMaxLR'] = impulsSHMMeanMaxLR
            dataByStim.loc[renderNames[ii], 'ImpulsSHM05ExMaxLR'] = impulsSHM05ExMaxLR

    # end of if branch for .asc files

    elif file[-4:] == "xlsx":
        workbookdata = pd.read_excel(io=file, sheet_name=None, engine='calamine')

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
        specFluctOldSHML = pd.DataFrame(workbookdata['Sheet6'].iloc[14:, 1:].values,
                                        columns=workbookdata['Sheet6'].iloc[13, 1:],
                                        index=workbookdata['Sheet6'].iloc[14:, 0])
        # right channel
        specFluctOldSHMR = pd.DataFrame(workbookdata['Sheet7'].iloc[14:, 1:].values,
                                        columns=workbookdata['Sheet7'].iloc[13, 1:],
                                        index=workbookdata['Sheet7'].iloc[14:, 0])

        # binaural specific fluctuation strength
        # (using ECMA-418-2:2022 Equation 112 for roughness)
        specFluctOldSHMBin = ((specFluctOldSHML**2
                               + specFluctOldSHMR**2)/2).pow(0.5)
        # binaural time-dependent fluctuation strength
        # (using ECMA-418-2:2022 Equation 111 for roughness)
        FluctOldSHMTDepBin = specFluctOldSHMBin.sum(axis=0)*bandDiff0p5

        # mask for start/end skip
        FluctOldSHMTDepBinMask = FluctOldSHMTDepBin.loc[(FluctOldSHMTDepBin.index.values
                                                         > start_skipT)
                                                        & (FluctOldSHMTDepBin.index.values
                                                           < FluctOldSHMTDepBin.index.values.max()
                                                           - end_skipT)]

        # binaural overall (90th percentile = 10% exceeded) fluctuation strength
        # (using ECMA-418-2:2022 Section 7.1.8 for roughness)
        FluctOldSHM10ExBin = FluctOldSHMTDepBinMask.quantile(q=0.90)

        # binaural overall (95th percentile = 5% exceeded) fluctuation strength
        FluctOldSHM05ExBin = FluctOldSHMTDepBinMask.quantile(q=0.95)

        # binaural overall time-averaged fluctuation strength
        FluctOldSHMAvgBin = FluctOldSHMTDepBinMask.mean()

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

        # max l/r overall (5% exceeded = 95th percentile) loudness
        loudISO105ExMaxLR = loudISO1TDepMask.quantile(q=0.95).max()

        # max l/r overall (power-averaged) loudness
        loudISO1PowAvgMaxLR = (loudISO1TDepMask.pow(1/np.log10(2)).mean(axis=0)**np.log10(2)).max()

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
        loudISO3PowAvgBin = (loudISO3TDepBinMask.pow(1/np.log10(2)).mean(axis=0)**np.log10(2)).iloc[0]

        # binaural overall 5% exceeded loudness
        loudISO305ExBin = loudISO3TDepBinMask.quantile(q=0.95).iloc[0]

        # Calculate overall DIN 45692 sharpness from 2-channel time-dependent
        # sharpness
        sharpDINTDep = pd.DataFrame(workbookdata['Sheet3'].iloc[13:, 1:3].values,
                                    columns=workbookdata['Sheet3'].iloc[12, 1:3],
                                    index=workbookdata['Sheet3'].iloc[13:, 0])

        # binaural sharpness
        sharpDINTDepBin = ((sharpDINTDep.iloc[:, 0]**2
                            + sharpDINTDep.iloc[:, 1]**2)/2).pow(0.5)

        # mask for start/end skip
        sharpDINTDepBinMask = sharpDINTDepBin.loc[(sharpDINTDepBin.index.values
                                                   > start_skipT)
                                                  & (sharpDINTDepBin.index.values
                                                     < sharpDINTDepBin.index.values.max()
                                                     - end_skipT)]

        # binaural overall (power-averaged) sharpness
        sharpDINPowAvgBin = sharpDINTDepBinMask.pow(1/np.log10(2)).mean()**np.log10(2)

        # binaural overall 5% exceeded sharpness
        sharpDIN05ExBin = sharpDINTDepBinMask.quantile(q=0.95)

        # Calculate overall von Bismarck | ISO 532-1 sharpness from 2-channel
        # time-dependent sharpness
        sharpvBISO1TDep = pd.DataFrame(workbookdata['Sheet4'].iloc[13:, 1:3].values,
                                       columns=workbookdata['Sheet4'].iloc[12, 1:3],
                                       index=workbookdata['Sheet4'].iloc[13:, 0])
        
        # binaural sharpness
        sharpvBISO1TDepBin = ((sharpvBISO1TDep.iloc[:, 0]**2
                               + sharpvBISO1TDep.iloc[:, 1]**2)/2).pow(0.5)

        # mask for start/end skip
        sharpvBISO1TDepBinMask = sharpvBISO1TDepBin.loc[(sharpvBISO1TDepBin.index.values
                                                         > start_skipT)
                                                        & (sharpvBISO1TDepBin.index.values
                                                           < sharpvBISO1TDepBin.index.values.max()    
                                                           - end_skipT)]

        # binaural overall (power-averaged) sharpness
        sharpvBISO1PowAvgBin = sharpvBISO1TDepBinMask.pow(1/np.log10(2)).mean()**np.log10(2)

        # binaural overall 5% exceeded sharpness
        sharpvBISO105ExBin = sharpvBISO1TDepBinMask.quantile(q=0.95)

        # Calculate overall Aures+ISO532-3 sharpness from 2-channel time-dependent
        # sharpness
        sharpAISO3TDep = pd.DataFrame(workbookdata['Sheet8'].iloc[13:, 1:3].values,
                                      columns=workbookdata['Sheet8'].iloc[12, 1:3],
                                      index=workbookdata['Sheet8'].iloc[13:, 0])

        # binaural sharpness
        sharpAISO3TDepBin = ((sharpAISO3TDep.iloc[:, 0]**2
                              + sharpAISO3TDep.iloc[:, 1]**2)/2).pow(0.5)

        # mask for start/end skip
        sharpAISO3TDepBinMask = sharpAISO3TDepBin.loc[(sharpAISO3TDepBin.index.values
                                                       > start_skipT)
                                                      & (sharpAISO3TDepBin.index.values
                                                         < sharpAISO3TDepBin.index.values.max()
                                                         - end_skipT)]

        # binaural overall (power-averaged) sharpness
        sharpAISO3PowAvgBin = sharpAISO3TDepBinMask.pow(1/np.log10(2)).mean()**np.log10(2)

        # binaural overall 5% exceeded sharpness
        sharpAISO305ExBin = sharpAISO3TDepBinMask.quantile(q=0.95)

        # Calculate overall Aures+ISO532-1 sharpness from 2-channel time-dependent
        # sharpness
        sharpAISO1TDep = pd.DataFrame(workbookdata['Sheet5'].iloc[13:, 1:3].values,
                                      columns=workbookdata['Sheet5'].iloc[12, 1:3],
                                      index=workbookdata['Sheet5'].iloc[13:, 0])

        # binaural sharpness
        sharpAISO1TDepBin = ((sharpAISO1TDep.iloc[:, 0]**2
                              + sharpAISO1TDep.iloc[:, 1]**2)/2).pow(0.5)

        # mask for start/end skip
        sharpAISO1TDepBinMask = sharpAISO1TDepBin.loc[(sharpAISO1TDepBin.index.values
                                                       > start_skipT)
                                                      & (sharpAISO1TDepBin.index.values
                                                         < sharpAISO1TDepBin.index.values.max()
                                                         - end_skipT)]

        # binaural overall (power-averaged) sharpness
        sharpAISO1PowAvgBin = sharpAISO1TDepBinMask.pow(1/np.log10(2)).mean()**np.log10(2)

        # binaural overall 5% exceeded sharpness
        sharpAISO105ExBin = sharpAISO1TDepBinMask.quantile(q=0.95)

        # binaural overall median sharpness
        sharpAISO1MedBin = sharpAISO1TDepBinMask.median()

        # add results to output DataFrame
        dataByStim.loc[renderNames[ii], 'FluctOldSHM10ExBin'] = FluctOldSHM10ExBin
        dataByStim.loc[renderNames[ii], 'FluctOldSHM05ExBin'] = FluctOldSHM05ExBin
        dataByStim.loc[renderNames[ii], 'FluctOldSHMAvgBin'] = FluctOldSHMAvgBin
        dataByStim.loc[renderNames[ii], 'LoudISO105ExMaxLR'] = loudISO105ExMaxLR
        dataByStim.loc[renderNames[ii], 'LoudISO1PowAvgMaxLR'] = loudISO1PowAvgMaxLR
        dataByStim.loc[renderNames[ii], 'LoudISO3PowAvgBin'] = loudISO3PowAvgBin
        dataByStim.loc[renderNames[ii], 'LoudISO305ExBin'] = loudISO305ExBin
        dataByStim.loc[renderNames[ii], 'SharpAurISO3PowAvgBin'] = sharpAISO3PowAvgBin
        dataByStim.loc[renderNames[ii], 'SharpAurISO305ExBin'] = sharpAISO305ExBin
        dataByStim.loc[renderNames[ii], 'SharpAurISO1PowAvgBin'] = sharpAISO1PowAvgBin
        dataByStim.loc[renderNames[ii], 'SharpAurISO105ExBin'] = sharpAISO105ExBin
        dataByStim.loc[renderNames[ii], 'SharpvBISO1PowAvgBin'] = sharpvBISO1PowAvgBin
        dataByStim.loc[renderNames[ii], 'SharpvBISO105ExBin'] = sharpvBISO105ExBin
        dataByStim.loc[renderNames[ii], 'SharpDINPowAvgBin'] = sharpDINPowAvgBin
        dataByStim.loc[renderNames[ii], 'SharpDIN05ExBin'] = sharpDIN05ExBin
        dataByStim.loc[renderNames[ii], 'SharpAISO1MedBin'] = sharpAISO1MedBin

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
            ambSharpAISO3TDepMovAvg = sharpAISO3TDep.rolling(window=int(np.ceil(sampleRateLoudISO3*windowT))).mean()

        elif renderNames[ii][0:3] in ["A1_", "A2_", "B2_"]:
            # calculate moving average values for combined stimulus
            sharpAISO3TDepMovAvg = sharpAISO3TDep.rolling(window=int(np.ceil(sampleRateLoudISO3*windowT))).mean()

            # calculate differences and make negative values 0
            dSharpAISO3TDep = np.maximum(sharpAISO3TDepMovAvg
                                         - ambSharpAISO3TDepMovAvg, 0)

            # binaural sharpness difference
            dSharpAISO3TDepBin = ((dSharpAISO3TDep.iloc[:, 0]**2
                                   + dSharpAISO3TDep.iloc[:, 1]**2)/2).pow(0.5)   

            # calculate aggregated difference values
            # binaural sharpness masked for start/end skip
            dSharpAISO3TDepBinMask = dSharpAISO3TDepBin.loc[(dSharpAISO3TDepBin.index.values
                                                             > start_skipT).transpose()
                                                            & (dSharpAISO3TDepBin.index.values
                                                               < dSharpAISO3TDepBin.index.values.max()
                                                               - end_skipT).transpose()]

            # binaural overall (power-averaged) sharpness
            dSharpAISO3PowAvgBin = dSharpAISO3TDepBinMask.pow(1/np.log10(2)).mean()**np.log10(2)

            # 2-channel overall 5% exceeded sharpness
            dSharpAISO305ExBin = dSharpAISO3TDepBinMask.quantile(q=0.95)

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'dSharpAurISO3PowAvgBin'] = dSharpAISO3PowAvgBin
            dataByStim.loc[renderNames[ii], 'dSharpAurISO305ExBin'] = dSharpAISO305ExBin

# end of for loop over ArtemiS SQM files

# --------------------------------------------------
# %% Psychoacoustic annoyance metrics single values
# --------------------------------------------------

dataByStim['PsychAnnoyWidmann'] = psych_annoy.widmannPA(dataByStim.loc[:, 'LoudISO105ExMaxLR'],
                                                        dataByStim.loc[:, 'SharpDIN05ExBin'],
                                                        dataByStim.loc[:, 'RoughFZ05ExMaxLR'],
                                                        dataByStim.loc[:, 'FluctFZ05ExMaxLR'])

dataByStim['PsychAnnoyMore'] = psych_annoy.morePA(dataByStim.loc[:, 'LoudISO105ExMaxLR'],
                                                  dataByStim.loc[:, 'SharpvBISO105ExBin'],
                                                  dataByStim.loc[:, 'RoughFZ05ExMaxLR'],
                                                  dataByStim.loc[:, 'FluctFZ05ExMaxLR'],
                                                  dataByStim.loc[:, 'TonalAur05ExMaxLR'])

dataByStim['PsychAnnoyDi'] = psych_annoy.diPA(dataByStim.loc[:, 'LoudISO105ExMaxLR'],
                                              dataByStim.loc[:, 'SharpvBISO105ExBin'],
                                              dataByStim.loc[:, 'RoughECMAAvgBin'],
                                              dataByStim.loc[:, 'FluctOldSHMAvgBin'],
                                              dataByStim.loc[:, 'TonalAurAvgMaxLR'])

dataByStim['PsychAnnoyTorija'] = psych_annoy.torijaPA(dataByStim.loc[:, 'LoudISO105ExMaxLR'],
                                                      dataByStim.loc[:, 'SharpDIN05ExBin'],
                                                      dataByStim.loc[:, 'RoughECMA05ExBin'],
                                                      dataByStim.loc[:, 'FluctECMA05ExBin'],
                                                      dataByStim.loc[:, 'TonalECMA05ExMaxLR'],
                                                      dataByStim.loc[:, 'ImpulsSHM05ExMaxLR'])

dataByStim['PsychAnnoyWillemsen'] = psych_annoy.willemsenPA(dataByStim.loc[:, 'LoudISO105ExMaxLR'],
                                                            dataByStim.loc[:, 'SharpAISO1MedBin'],
                                                            dataByStim.loc[:, 'RoughECMA05ExBin'],
                                                            dataByStim.loc[:, 'ImpulsLoudWZAvgMaxLR'])

dataByStim['PsychAnnoyBoucher'] = psych_annoy.boucherPA(dataByStim.loc[:, 'LoudISO105ExMaxLR'],
                                                        dataByStim.loc[:, 'SharpDIN05ExBin'],
                                                        dataByStim.loc[:, 'RoughECMAAvgBin'],
                                                        dataByStim.loc[:, 'FluctOldSHMAvgBin'],
                                                        dataByStim.loc[:, 'TonalECMAAvgMaxLR'])


# for rows in dataByStim.index that start with "A1_", (but not "A1"), fill in the dPsychAnnoy... columns by
# subtracting the corresponding "A1" row values from the "A1_" row values 
for row in dataByStim.index:
    if row.startswith("A1_") and row != "A1":
        base_row = "A1"
        for col in ['dPsychAnnoyWidmann', 'dPsychAnnoyMore', 'dPsychAnnoyDi', 'dPsychAnnoyTorija', 'dPsychAnnoyWillemsen', 'dPsychAnnoyBoucher']:
            dataByStim.loc[row, col] = dataByStim.loc[row, col.replace('dPsychAnnoy', 'PsychAnnoy')] - dataByStim.loc[base_row, col.replace('dPsychAnnoy', 'PsychAnnoy')]
    elif row.startswith("A2_") and row != "A2":
        base_row = "A2"
        for col in ['dPsychAnnoyWidmann', 'dPsychAnnoyMore', 'dPsychAnnoyDi', 'dPsychAnnoyTorija', 'dPsychAnnoyWillemsen', 'dPsychAnnoyBoucher']:
            dataByStim.loc[row, col] = dataByStim.loc[row, col.replace('dPsychAnnoy', 'PsychAnnoy')] - dataByStim.loc[base_row, col.replace('dPsychAnnoy', 'PsychAnnoy')]
    elif row.startswith("B2_") and row != "B2":
        base_row = "B2"
        for col in ['dPsychAnnoyWidmann', 'dPsychAnnoyMore', 'dPsychAnnoyDi', 'dPsychAnnoyTorija', 'dPsychAnnoyWillemsen', 'dPsychAnnoyBoucher']:
            dataByStim.loc[row, col] = dataByStim.loc[row, col.replace('dPsychAnnoy', 'PsychAnnoy')] - dataByStim.loc[base_row, col.replace('dPsychAnnoy', 'PsychAnnoy')]


# add UAS-only and ambient-only data to combined stimuli
# ------------------------------------------------------
indicesDiffPsycho = ["PartLoudSHMPowAvg",
                     "PartTonLdSHMPowAvg",
                     "PartTonShpAurSHMPowAvg",
                     "PartTonShpAurSHM05Ex",
                     "dTonalECMAAvg",
                     "dTonalECMA05Ex",
                     "dTonalSHMIntAvg",
                     "dTonalSHMInt05Ex",
                     "dTonalAwSHMAvg",
                     "dTonalAwSHM05Ex",
                     "dTonalAwSHMIntAvg",
                     "dTonalAwSHMInt05Ex",
                     "dTonLdECMAPowAvg",
                     "dTonLdECMA05Ex",
                     "dTonShpAurSHMPowAvg",
                     "dTonShpAurSHM05Ex",
                     "dRoughECMA10Ex",
                     "dRoughECMA05Ex",
                     "dRoughFZ10Ex",
                     "dRoughFZ05Ex",
                     "dFluctECMA10Ex",
                     "dFluctECMA05Ex",
                     "dFluctOV10Ex",
                     "dFluctOV05Ex",
                     "dSharpAurSHMPowAvg",
                     "dSharpAurSHM05Ex",
                     "dSharpAurISO3PowAvg",
                     "dSharpAurISO305Ex",
                     "dImpulsSHMPowAvg",
                     "dImpulsSHMAvg",
                     "dImpulsSHM05Ex",
                     "dImpulsLoudWZAvg",
                     "dImpulsLoudWZ05Ex",
                     "dImpulsLoudWZPowAvg",
                     "dImpulsLoudWECMAAvg",
                     "dImpulsLoudWECMA05Ex",
                     "dImpulsLoudWECMAPowAvg",
                     "dPsychAnnoyWidmann",
                     "dPsychAnnoyMore",
                     "dPsychAnnoyDi",
                     "dPsychAnnoyTorija",
                     "dPsychAnnoyWillemsen",
                     "dPsychAnnoyBoucher"]

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

indicesLevelDiffs = ['LAeqLAF90diff', 'LASmaxLAF90diff', 'LASmaxLAF50diff',
                     'LASmaxLAeqdiff', 'LAELAF90diff', 'LAELAF50diff',
                     'LAELAeqdiff', 'PNLMLAeqdiff', 'PNLTMLAeqdiff',
                     'EPNLLAeqdiff', 'PNLMLAF90diff', 'PNLTMLAF90diff',
                     'EPNLLAF90diff', 'PNLMLAF50diff', 'PNLTMLAF50diff',
                     'EPNLLAF50diff']

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
                            'PartLoudMGLTPowAvg', 'PartLoudMGLT05Ex',
                            'PartLoudSHMPowAvgBin',] ] = 0

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
# -----------------------

# ensure level difference metrics follow absolute SQMs
indicesLevelDiffs.reverse()
cols = list(dataByStim.columns)
for col in indicesLevelDiffs:
    cols.insert(cols.index(indicesAbsPsycho[-1])
                + 1, cols.pop(cols.index(col)))
dataByStim = dataByStim[cols]
indicesLevelDiffs.reverse()  # revert list for use later

# move partial tonal metrics to after partial loudness metrics
cols = list(dataByStim.columns)
cols.insert(cols.index(partLoudness.columns[-1]) + 1,
            cols.pop(cols.index('PartLoudSHMPowAvgBin')))
cols.insert(cols.index('PartLoudSHMPowAvgBin') + 1,
            cols.pop(cols.index('PartTonLdSHMPowAvgBin')))
cols.insert(cols.index('PartTonLdSHMPowAvgBin') + 1,
            cols.pop(cols.index('PartTonShpAurSHMPowAvgBin')))
cols.insert(cols.index('PartTonShpAurSHMPowAvgBin') + 1,
            cols.pop(cols.index('PartTonShpAurSHM05ExBin')))
dataByStim = dataByStim[cols]

# move detectability metrics to after partial loudness metrics
indicesDetect.reverse()
cols = list(dataByStim.columns)
for col in indicesDetect:
    cols.insert(cols.index(partLoudness.columns[-1]) + 1,
                cols.pop(cols.index(col)))
dataByStim = dataByStim[cols]
indicesDetect.reverse()  # revert list for use later   

# move SQM difference metrics to after detectability metrics
indicesDiffPsycho.reverse()
cols = list(dataByStim.columns)
for col in indicesDiffPsycho:
    cols.insert(cols.index(indicesDetect[-1]) + 1,
                cols.pop(cols.index(col)))
dataByStim = dataByStim[cols]
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
                               header=0, engine='calamine')
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

# loop over stimuli recording names, extract corresponding response data and
# tranpose data with recording as index
# calculate aggregate statistics
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
    else:
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
                               header=0, engine='calamine')
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
    else:
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

# merge response data into output
dataByStim = dataByStim.merge(allResponses, how='outer',
                              left_index=True, right_index=True)

# add pre- and post-test response responses to dataset
preTestResponses = pd.read_excel(io=filepath, sheet_name="Pre-test",
                                 header=0, engine='calamine')
preTestResponses.drop(columns=['Impairment_details'], inplace=True)
preTestResponses.drop(columns=preTestResponses.loc[:,
                                                   'PANAS_interested':].columns,
                      inplace=True)
postTestResponses = pd.read_excel(io=filepath, sheet_name="Post-test",
                                  header=0, engine='calamine')
postTestResponses.drop(columns=['Nationality', 'Language',
                                'Response Considerations',
                                'Source comments', 'AAM Exp',
                                'General Feedback'], inplace=True)
postTestResponses.drop(columns=postTestResponses.loc[:, 'NSS21':'NSS8'].columns,
                       inplace=True)
# replace problem columns names
postTestResponses.columns = postTestResponses.columns.str.replace(" ", "_")

prePostTestResponses = preTestResponses.merge(postTestResponses, on='ID#')

# set Exp2ID column to pd.int 64 type
prePostTestResponses['Exp2ID'] = prePostTestResponses['Exp2ID'].astype(pd.Int64Dtype())

# ----------------------------------
# Prepare outputs for saving to file
# ----------------------------------

# separate 'by stimulus' output into test data and auxiliary data and save to
# file
dataByStimTest = dataByStim.loc[(dataByStim.index.str.find("B2") == 0)
                                | (dataByStim.index.str.find("A1") == 0)
                                | (dataByStim.index.str.find("A2") == 0), :]
# set stimID to integer type
dataByStimTest['StimID'] = dataByStimTest['StimID'].astype(int)

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

# merge response and stimuli data into 'by participant' test datasets, and
# save to file

partADataBySubj = pd.merge(left=partAData,
                           right=dataByStimTestA.loc[:, :dataByStimTestA.columns[dataByStimTestA.columns.get_loc('Valence_1') - 1]],
                           how='outer', left_on='Recording', right_on='CALBINRecFiles')
partADataBySubj.sort_values(by='ID#', axis=0, inplace=True)
partADataBySubj = pd.merge(left=partADataBySubj,
                           right=prePostTestResponses, how='left',
                           left_on='ID#', right_on='ID#')
partADataBySubj.drop(columns=['Part', 'Recording'], inplace=True)
partADataBySubj.insert(loc=0, column='SessionPart',
                       value=partADataBySubj.pop('SessionPart'))

partBDataBySubj = pd.merge(left=partBData,
                           right=dataByStimTestB.loc[:, :dataByStimTestB.columns[dataByStimTestB.columns.get_loc('Valence_1') - 1]],
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