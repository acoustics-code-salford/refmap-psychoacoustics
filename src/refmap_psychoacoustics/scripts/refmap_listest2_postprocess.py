# -*- coding: utf-8 -*-

# script


# %%%%%
# Setup
# -----

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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Initialise data frame for metric calculations and results
# ---------------------------------------------------------

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
filelist = list(QFileDialog.getOpenFileNames(caption="Open recording files in '03 Experiment\Experiment 2\Stimuli\Calibrated_recordings\RecordHATS'",
                                             filter=fileExts))[0]
filelist.sort()
filenames = [filepath.split('/')[-1] for filepath in filelist]
stemNames = [filename.replace("_HATS_Pa.wav", "") for filename in filenames]
sqmNames = [filename.replace("HATS", "MA220") for filename in filenames
            if filename.find("G57") == -1]


dataByStim = pd.DataFrame(index=stemNames, dtype=float)

dataByStim['HATSRecFiles'] = filenames
dataByStim['MA220MicRecFiles'] = sqmNames

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Add categorical variables for stimuli
# -------------------------------------

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

# ambient env ref
dataByStim.loc[dataByStim.index.str.find("Background") == -1, 'AmbientRef'] = "UAS only"
dataByStim.loc[dataByStim.index.str.find("BusyStreet6") != -1, 'AmbientRef'] = "BusyStreet6"
dataByStim.loc[dataByStim.index.str.find("BusyStreet6Quiet") != -1, 'AmbientRef'] = "BusyStreet6Quiet"
dataByStim.loc[dataByStim.index.str.find("BusyStreet8") != -1, 'AmbientRef'] = "BusyStreet8"
dataByStim.loc[dataByStim.index.str.find("Park3") != -1, 'AmbientRef'] = "Park3"
dataByStim.loc[dataByStim.index.str.find("Park3Loud") != -1, 'AmbientRef'] = "Park3Loud"
dataByStim.loc[dataByStim.index.str.find("QuietStreet7") != -1, 'AmbientRef'] = "QuietStreet7"

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

# UAS proximity
dataByStim.loc[dataByStim.index.str.find("Near") != -1, 'UASProximity'] = "Near"
dataByStim.loc[dataByStim.index.str.find("Far") != -1, 'UASProximity'] = "Far"
dataByStim.loc[(dataByStim.index.str.find("Baseline") != -1), 'UASProximity'] = "Baseline"

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
dataByStim.loc[dataByStim.index.str.find("Baseline") != -1, 'UASStart'] = "Baseline"

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    signalA = filterFuncs.A_weight_T(signal, sampleRatein)
    signalmagAF = filterFuncs.time_weight(signalA, sampleRatein, tau=0.125)
    signalmagAS = filterFuncs.time_weight(signalA, sampleRatein, tau=1)

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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Psychoacoustic metrics single values calculation
# ------------------------------------------------

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
                 "PartSharpAurSHMPowAvg",
                 "PartSharpAurSHM05Ex",
                 "PartSharpWidSHMPowAvg",
                 "PartSharpWidSHM05Ex",
                 "PartSharpvBSHMPowAvg",
                 "PartSharpvBSHM05Ex",
                 "PartTonShpAurSHMPowAvg",
                 "PartTonShpAurSHM05Ex",
                 "PartTonShpWidSHMPowAvg",
                 "PartTonShpWidSHM05Ex",
                 "PartTonShpvBSHMPowAvg",
                 "PartTonShpvBSHM05Ex",
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

# %%%%%%%%%%%%%%%%%%%%%%%%
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
renderNames = [filename.split('_MA220_Pa')[0] for filename in filenames]

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
            fluctFZ05Ex = fluctFZTDepMask.quantile(q=0.95).iloc[0]

            # overall 10% exceeded fluctuation strength
            fluctFZ10Ex = fluctFZTDepMask.quantile(q=0.90).iloc[0]
            
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
            roughFZ05Ex = roughFZTDepMask.quantile(q=0.95).iloc[0]

            # overall 10% exceeded roughness
            roughFZ10Ex = roughFZTDepMask.quantile(q=0.90).iloc[0]

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
            tonLdECMA05Ex = tonLdECMATDepMask.quantile(q=0.95)
            
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
            tonalECMAAvg = tonalECMATDepMask.mean()

            # 5% exceeded tonality (omitting T<=0.02)
            tonalECMA05Ex = tonalECMATDepMask.quantile(q=0.95)

            # time-averaged integrated tonality (omitting T<=0.02)
            # NOTE: uses mask from ECMA tonality
            tonalSHMIntAvg = tonalSHMIntTDepMask.mean()

            # 5% exceeded integrated tonality (omitting T<=0.02)
            tonalSHMInt05Ex = tonalSHMIntTDepMask.quantile(q=0.95)

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
            tonalAwSHMAvg = tonalAwSHMTDepMask.mean()

            # 5% exceeded tonal annoyance-weighted tonality (omitting T<=0.02)
            tonalAwSHM05Ex = tonalAwSHMTDepMask.quantile(q=0.95)

            # time-averaged integrated tonal annoyance-weighted tonality (omitting T<=0.02)
            # NOTE: uses mask from ECMA tonal annoyance-weighted tonality
            tonalAwSHMIntAvg = tonalAwSHMIntTDepMask.mean()

            # 5% exceeded integrated tonal annoyance-weighted tonality (omitting T<=0.02)
            tonalAwSHMInt05Ex = tonalAwSHMIntTDepMask.quantile(q=0.95)

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
            roughECMA10Ex = roughECMATDepMask.quantile(q=0.90)

            # overall (95th percentile = 5% exceeded) roughness
            roughECMA05Ex = roughECMATDepMask.quantile(q=0.95)

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
            fluctOV10Ex = fluctOVTDepMask.quantile(q=0.90)

            # overall (95th percentile = 5% exceeded) fluctuation strength
            fluctOV05Ex = fluctOVTDepMask.quantile(q=0.95)

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
                                                           - end_skipT)]
            
            # overall (averaged) impulsive loudness
            impulsLoudWZAvg = impulsLoudWZTDepMask.mean().iloc[0]

            # overall (95th percentile = 5% exceeded) impulsive loudness
            impulsLoudWZ05Ex = impulsLoudWZTDepMask.quantile(q=0.95).iloc[0]

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
                                                                 - end_skipT)]

            # overall (averaged) impulsive loudness
            impulsLoudWECMAAvg = impulsLoudWECMATDepMask.mean().iloc[0]

            # overall (95th percentile = 5% exceeded) impulsive loudness
            impulsLoudWECMA05Ex = impulsLoudWECMATDepMask.quantile(q=0.95).iloc[0]

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
            partLoudSHMPowAvg = partLoudSHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)

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
            partTonLdSHMPowAvg = partTonLdSHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)

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
            if renderNames[ii].__contains__("_Baseline"):
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
                    dRoughFZ10Ex = dRoughFZTDepMask.quantile(q=0.90).iloc[0]

                    # overall (95th percentile = 5% exceeded) roughness
                    dRoughFZ05Ex = dRoughFZTDepMask.quantile(q=0.95).iloc[0]

                    # add results to output DataFrame
                    dataByStim.loc[renderNames[ii], 'dRoughFZ10Ex'] = dRoughFZ10Ex
                    dataByStim.loc[renderNames[ii], 'dRoughFZ05Ex'] = dRoughFZ05Ex

                elif filenames[ii].find("TonalECMA") != -1:
                    specTonalECMAMovAvg = specTonalECMA.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean()
                    
                    # calculate differences and make negative values 0
                    dSpecTonalECMA = np.maximum(specTonalECMAMovAvg
                                                - ambSpecTonalECMAMovAvg, 0)

                    # time-dependent tonality (max, not integration)
                    dTonalECMATDep = dSpecTonalECMA.max(axis=1).to_frame(name="Monaural")
                    
                    # time-dependent integrated tonality
                    dTonalSHMIntTDep = (dSpecTonalECMA.sum(axis=1)*bandDiff0p5*0.348088948583815).to_frame(name="Monaural")

                    # mask for start/end skip and values <= 0.02
                    dTonalECMATDepMask = dTonalECMATDep.loc[(dTonalECMATDep.index.values
                                                             > start_skipT)
                                                            & (dTonalECMATDep.index.values
                                                               < dTonalECMATDep.index.values.max()
                                                               - end_skipT)
                                                            & (dTonalECMATDep.values
                                                               > 0.02)]

                    # mask for start/end skip and values <= 0.02
                    # NOTE: uses tonality magnitude mask from ECMA tonality
                    dTonalSHMIntTDepMask = dTonalSHMIntTDep.loc[(dTonalSHMIntTDep.index.values
                                                                 > start_skipT)
                                                                & (dTonalSHMIntTDep.index.values
                                                                   < dTonalSHMIntTDep.index.values.max()
                                                                   - end_skipT)
                                                                & (dTonalECMATDep.values
                                                                   > 0.02)]  # see NOTE above

                    # time-averaged tonality (omitting T<=0.02)
                    dTonalECMAAvg = dTonalECMATDepMask.mean().iloc[0]

                    # 5% exceeded tonality (automatically omitting T<=0.02)
                    dTonalECMA05Ex = dTonalECMATDepMask.quantile(q=0.95).iloc[0]

                    # time-averaged integrated tonality (omitting T<=0.02)
                    dTonalSHMIntAvg = dTonalSHMIntTDepMask.mean().iloc[0]
                    # 5% exceeded integated tonality (automatically omitting T<=0.02)
                    dTonalSHMInt05Ex = dTonalSHMIntTDepMask.quantile(q=0.95).iloc[0]

                    # add results to output DataFrame
                    dataByStim.loc[renderNames[ii], 'dTonalECMAAvg'] = dTonalECMAAvg
                    dataByStim.loc[renderNames[ii], 'dTonalECMA05Ex'] = dTonalECMA05Ex
                    dataByStim.loc[renderNames[ii], 'dTonalSHMIntAvg'] = dTonalSHMIntAvg
                    dataByStim.loc[renderNames[ii], 'dTonalSHMInt05Ex'] = dTonalSHMInt05Ex

                    # tonal annoyance-weighted
                    specTonalAwSHMMovAvg = specTonalAwSHM.rolling(window=int(np.ceil(sampleRateTonalECMA*windowT))).mean()
        
                    # calculate differences and make negative values 0
                    dSpecTonalAwSHM = np.maximum(specTonalAwSHMMovAvg
                                                 - ambSpecTonalAwSHMMovAvg, 0)

                    # time-dependent tonal annoyance-weighted tonality
                    dTonalAwSHMTDep = dSpecTonalAwSHM.max(axis=1).to_frame(name="Monaural")

                    # time-dependent integrated tonal annoyance-weighted tonality
                    dTonalAwSHMIntTDep = (dSpecTonalAwSHM.sum(axis=1)*bandDiff0p5*0.348088948583815).to_frame(name="Monaural")
                    
                    # mask for start/end skip and values <= 0.02
                    dTonalAwSHMTDepMask = dTonalAwSHMTDep.loc[(dTonalAwSHMTDep.index.values
                                                               > start_skipT)
                                                              & (dTonalAwSHMTDep.index.values
                                                                 < dTonalAwSHMTDep.index.values.max()
                                                                 - end_skipT)
                                                              & (dTonalAwSHMTDep.values
                                                                 > 0.02)]

                    # mask for start/end skip and values <= 0.02
                    # NOTE: uses tonality magnitude mask from ECMA tonality
                    dTonalAwSHMIntTDepMask = dTonalAwSHMIntTDep.loc[(dTonalAwSHMIntTDep.index.values
                                                                     > start_skipT)
                                                                    & (dTonalAwSHMIntTDep.index.values
                                                                       < dTonalAwSHMIntTDep.index.values.max()
                                                                       - end_skipT)
                                                                    & (dTonalAwSHMTDep.values
                                                                       > 0.02)]  # see NOTE above

                    # time-averaged tonality (omitting T<=0.02)
                    dTonalAwSHMAvg = dTonalAwSHMTDepMask.mean().iloc[0]

                    # 5% exceeded tonality (automatically omitting T<=0.02)
                    dTonalAwSHM05Ex = dTonalAwSHMTDepMask.quantile(q=0.95).iloc[0]

                    # time-averaged integrated tonality (omitting T<=0.02)
                    dTonalAwSHMIntAvg = dTonalAwSHMIntTDepMask.mean().iloc[0]

                    # 5% exceeded integrated tonality (automatically omitting T<=0.02)
                    dTonalAwSHMInt05Ex = dTonalAwSHMIntTDepMask.quantile(q=0.95).iloc[0]

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
                    dTonLdECMA05Ex = dTonLdECMATDepMask.quantile(q=0.95)

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
                    dRoughECMA10Ex = dRoughECMATDepMask.quantile(q=0.90)
                    
                    # overall (95th percentile = 5% exceeded) roughness
                    dRoughECMA05Ex = dRoughECMATDepMask.quantile(q=0.95)

                    # add results to output DataFrame
                    dataByStim.loc[renderNames[ii], 'dRoughECMA10Ex'] = dRoughECMA10Ex
                    dataByStim.loc[renderNames[ii], 'dRoughECMA05Ex'] = dRoughECMA05Ex

                elif filenames[ii].find("FluctOV") != -1:
                    specFluctOVMovAvg = specFluctOV.rolling(window=int(np.ceil(sampleRateFluctOV*windowT))).mean()

                    # calculate differences and make negative values 0
                    dSpecFluctOV = np.maximum(specFluctOVMovAvg
                                              - ambSpecFluctOVMovAvg, 0)

                    # time-dependent integrated fluctuation strength
                    dFluctOVTDep = dSpecFluctOV.sum(axis=1)*bandDiff0p5
                    
                    # fluctuation strength masked for start/end skip
                    dFluctOVTDepMask = dFluctOVTDep.loc[(dFluctOVTDep.index.values
                                                         > start_skipT).transpose()
                                                        & (dFluctOVTDep.index.values
                                                           < dFluctOVTDep.index.values.max()
                                                           - end_skipT).transpose()]

                    # overall 5% exceeded fluctuation strength 
                    dFluctOV05Ex = dFluctOVTDepMask.quantile(q=0.95)

                    # overall 10% exceeded fluctuation strength 
                    dFluctOV10Ex = dFluctOVTDepMask.quantile(q=0.90)
                    
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
                                                                     - end_skipT)]

                    # overall (averaged) impulsive loudness
                    dImpulsLoudWZAvg = dImpulsLoudWZTDepMask.mean().iloc[0]
                    
                    # overall (95th percentile = 5% exceeded) impulsive loudness
                    dImpulsLoudWZ05Ex = dImpulsLoudWZTDepMask.quantile(q=0.95).iloc[0]
                    
                    # overall (power-averaged) impulsive loudness
                    dImpulsLoudWZPowAvg = (dImpulsLoudWZTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

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
                                                                           - end_skipT)]

                    # overall (averaged) impulsive loudness
                    dImpulsLoudWECMAAvg = dImpulsLoudWECMATDepMask.mean().iloc[0]
                    
                    # overall (95th percentile = 5% exceeded) impulsive loudness
                    dImpulsLoudWECMA05Ex = dImpulsLoudWECMATDepMask.quantile(q=0.95, axis=0).iloc[0]
                    
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
        loudQZ5321PowAvg = (loudQZ5321TDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

        # overall 5% exceeded loudness
        loudQZ532105Ex = loudQZ5321TDepMask.quantile(q=0.95).iloc[0]
        
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
        loudQZ5323PowAvg = (loudQZ5323TDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

        # overall 5% exceeded loudness
        loudQZ532305Ex = loudQZ5323TDepMask.quantile(q=0.95).iloc[0]

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
        loudQZ4182PowAvg = (loudQZ4182TDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]
    
        # overall 5% exceeded loudness
        loudQZ418205Ex = loudQZ4182TDepMask.quantile(q=0.95).iloc[0]

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
        sharpAQZ532105Ex = sharpAQZ5321TDepMask.quantile(q=0.95).iloc[0]

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
        sharpAQZ532305Ex = sharpAQZ5323TDepMask.quantile(q=0.95).iloc[0]
        
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
        sharpAQZ418205Ex = sharpAQZ4182TDepMask.quantile(q=0.95).iloc[0]
        
        # Calculate overall Aures tonality from time-dependent
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
        tonalAurAvg = tonalAurTDepMask.mean().iloc[0]

        # overall 5% exceeded tonality
        tonalAur05Ex = tonalAurTDepMask.quantile(q=0.95).iloc[0]

        # overall 10% exceeded tonality
        tonalAur10Ex = tonalAurTDepMask.quantile(q=0.90).iloc[0]

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
        roughDW10Ex = roughDWTDepMask.quantile(q=0.90)

        # overall (95th percentile = 5% exceeded) roughness
        roughDW05Ex = roughDWTDepMask.quantile(q=0.95)

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
        sharpASHMPowAvg = (sharpASHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

        # overall 5% exceeded sharpness
        sharpASHM05Ex = sharpASHMTDepMask.quantile(q=0.95).iloc[0]

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
        tonShpASHM05Ex = tonShpASHMTDepMask.quantile(q=0.95).iloc[0]

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

        # calculation section for Aures+Sottek Hearing Model partial sharpness
        if workbookdata.keys().__contains__('PartShpAuresSHM'):

            # Calculate overall Aures+Sottek Hearing Model partial tonal sharpness from 2-channel
            # time-dependent partial sharpness
            partSharpASHMTDep = pd.DataFrame(workbookdata['PartShpAuresSHM'].iloc[0:, 1:].values,
                                             columns=workbookdata['PartShpAuresSHM'].iloc[0, 1:].index,
                                             index=workbookdata['PartShpAuresSHM'].iloc[:, 0])

            # mask for start/end skip
            partSharpASHMTDepMask = partSharpASHMTDep.loc[(partSharpASHMTDep.index.values
                                                           > start_skipT)
                                                          & (partSharpASHMTDep.index.values
                                                             < partSharpASHMTDep.index.values.max()
                                                             - end_skipT)]

            # overall (power-averaged) partial sharpness
            partSharpASHMPowAvg = (partSharpASHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall 5% exceeded partial sharpness
            partSharpASHM05Ex = partSharpASHMTDepMask.quantile(q=0.95).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'PartSharpAurSHMPowAvg'] = partSharpASHMPowAvg
            dataByStim.loc[renderNames[ii], 'PartSharpAurSHM05Ex'] = partSharpASHM05Ex
        
        # calculation section for Widmann+Sottek Hearing Model partial sharpness
        if workbookdata.keys().__contains__('PartShpDINSHM'):

            # Calculate overall Widmann+Sottek Hearing Model partial tonal sharpness from 2-channel
            # time-dependent partial sharpness
            partSharpWSHMTDep = pd.DataFrame(workbookdata['PartShpDINSHM'].iloc[0:, 1:].values,
                                             columns=workbookdata['PartShpDINSHM'].iloc[0, 1:].index,
                                             index=workbookdata['PartShpDINSHM'].iloc[:, 0])

            # mask for start/end skip
            partSharpWSHMTDepMask = partSharpWSHMTDep.loc[(partSharpWSHMTDep.index.values
                                                           > start_skipT)
                                                          & (partSharpWSHMTDep.index.values
                                                             < partSharpWSHMTDep.index.values.max()
                                                             - end_skipT)]

            # overall (power-averaged) partial sharpness
            partSharpWSHMPowAvg = (partSharpWSHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall 5% exceeded partial sharpness
            partSharpWSHM05Ex = partSharpWSHMTDepMask.quantile(q=0.95).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'PartSharpWidSHMPowAvg'] = partSharpWSHMPowAvg
            dataByStim.loc[renderNames[ii], 'PartSharpWidSHM05Ex'] = partSharpWSHM05Ex

        # calculation section for von Bismarck+Sottek Hearing Model partial sharpness
        if workbookdata.keys().__contains__('PartShpvBSHM'):

            # Calculate overall von Bismarck+Sottek Hearing Model partial tonal sharpness from 2-channel
            # time-dependent partial sharpness
            partSharpvBSHMTDep = pd.DataFrame(workbookdata['PartShpvBSHM'].iloc[0:, 1:].values,
                                              columns=workbookdata['PartShpvBSHM'].iloc[0, 1:].index,
                                              index=workbookdata['PartShpvBSHM'].iloc[:, 0])
            # mask for start/end skip
            partSharpvBSHMTDepMask = partSharpvBSHMTDep.loc[(partSharpvBSHMTDep.index.values
                                                             > start_skipT)
                                                            & (partSharpvBSHMTDep.index.values
                                                               < partSharpvBSHMTDep.index.values.max()
                                                               - end_skipT)]

            # overall (power-averaged) partial sharpness
            partSharpvBSHMPowAvg = (partSharpvBSHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall 5% exceeded partial sharpness
            partSharpvBSHM05Ex = partSharpvBSHMTDepMask.quantile(q=0.95).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'PartSharpvBSHMPowAvg'] = partSharpvBSHMPowAvg
            dataByStim.loc[renderNames[ii], 'PartSharpvBSHM05Ex'] = partSharpvBSHM05Ex

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
            partTonShpASHMPowAvg = (partTonShpASHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall 5% exceeded partial tonal sharpness
            partTonShpASHM05Ex = partTonShpASHMTDepMask.quantile(q=0.95).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'PartTonShpAurSHMPowAvg'] = partTonShpASHMPowAvg
            dataByStim.loc[renderNames[ii], 'PartTonShpAurSHM05Ex'] = partTonShpASHM05Ex
        
        # calculation section for Widmann+Sottek Hearing Model partial tonal sharpness        
        if workbookdata.keys().__contains__('PartTonShpDINSHM'):

            # Calculate overall Widmann+Sottek Hearing Model partial tonal sharpness from
            # time-dependent partial tonal sharpness
            partTonShpWSHMTDep = pd.DataFrame(workbookdata['PartTonShpDINSHM'].iloc[0:, 1:].values,
                                              columns=workbookdata['PartTonShpDINSHM'].iloc[0, 1:].index,
                                              index=workbookdata['PartTonShpDINSHM'].iloc[:, 0])

            # mask for start/end skip
            partTonShpWSHMTDepMask = partTonShpWSHMTDep.loc[(partTonShpWSHMTDep.index.values
                                                             > start_skipT)
                                                            & (partTonShpWSHMTDep.index.values
                                                               < partTonShpWSHMTDep.index.values.max()
                                                               - end_skipT)]

            # overall (power-averaged) partial tonal sharpness
            partTonShpWSHMPowAvg = (partTonShpWSHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall 5% exceeded partial tonal sharpness
            partTonShpWSHM05Ex = partTonShpWSHMTDepMask.quantile(q=0.95).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'PartTonShpWidSHMPowAvg'] = partTonShpWSHMPowAvg
            dataByStim.loc[renderNames[ii], 'PartTonShpWidSHM05Ex'] = partTonShpWSHM05Ex

        # calculation section for von Bismarck+Sottek Hearing Model partial tonal sharpness        
        if workbookdata.keys().__contains__('PartTonShpvBSHM'):

            # Calculate overall von Bismarck+Sottek Hearing Model partial tonal sharpness from
            # time-dependent partial tonal sharpness
            partTonShpvBSHMTDep = pd.DataFrame(workbookdata['PartTonShpvBSHM'].iloc[0:, 1:].values,
                                               columns=workbookdata['PartTonShpvBSHM'].iloc[0, 1:].index,
                                               index=workbookdata['PartTonShpvBSHM'].iloc[:, 0])

            # mask for start/end skip
            partTonShpvBSHMTDepMask = partTonShpvBSHMTDep.loc[(partTonShpvBSHMTDep.index.values
                                                               > start_skipT)
                                                              & (partTonShpvBSHMTDep.index.values
                                                                 < partTonShpvBSHMTDep.index.values.max()
                                                                 - end_skipT)]

            # overall (power-averaged) partial tonal sharpness
            partTonShpvBSHMPowAvg = (partTonShpvBSHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall 5% exceeded partial tonal sharpness
            partTonShpvBSHM05Ex = partTonShpvBSHMTDepMask.quantile(q=0.95).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'PartTonShpvBSHMPowAvg'] = partTonShpvBSHMPowAvg
            dataByStim.loc[renderNames[ii], 'PartTonShpvBSHM05Ex'] = partTonShpvBSHM05Ex

        # calculation section for SQM differences
        # NOTE: THIS SECTION RELIES ON THE ALPHABETIC ORDER OF THE STIMULI FILES AS
        # ORIGINALLY NAMED: EACH AMBIENT SOUND FILE PRECEDING THE CORRESPONDING
        # COMBINED STIMULI FILES. IF THE ORDERING OR FILE NAMING IS CHANGED, THIS
        # CALCULATION WILL BE INVALIDATED AND *MAY ALSO* CAUSE AN ERROR.
        # a rolling 50 ms window is applied to average the SQM values over time -
        # this is to reduce uncertainty due to imperfect time-alignment between the
        # ambient vs combined stimuli files (all are from recordings, so there will
        # be some slippage due to imperfect editing)
        if renderNames[ii].__contains__("_Baseline"):
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
            dSharpASHM05Ex = dSharpASHMTDepMask.quantile(q=0.95).iloc[0]

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
            dTonShpASHM05Ex = dTonShpASHMTDepMask.quantile(q=0.95).iloc[0]
            
            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'dTonShpAurSHMPowAvg'] = dTonShpASHMPowAvg
            dataByStim.loc[renderNames[ii], 'dTonShpAurSHM05Ex'] = dTonShpASHM05Ex

# end of for loop over MATLAB SQM files

# %%%%%%%%%%%%%%%%%%%%%%%%%
# From ArtemiS calculations
# -------------------------

# open xlsx file selection dialog and assign filepaths to list
# PROJECT NOTE: the calculated files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%202/Analysis/ArtemiS?csf=1&web=1&e=mIwerm
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance()

fileExts = "*.asc"
filelist = list(QFileDialog.getOpenFileNames(filter=fileExts,
                                             caption=r"Select ArtemiS output files in '03 Experiment\Experiment 2\Analysis\ArtemiS\asc'"))[0]

filelist.sort()

filenames = [filepath.split('/')[-1] for filepath in filelist]
renderNames = [filename.split('_MA220_Pa')[0] for filename in filenames]

# loop over files to analyse
for ii, file in enumerate(filelist):
    if ii == 0:
        print("Processing sound quality metrics...\n")
    print(file.split('/')[-1])

    if file.find("Specific") != -1:

        if file.find("Impulsiveness") != -1:
            # Sottek Hearing Model impulsiveness
            # ----------------------------------
            filedata = np.array(pd.read_csv(file, sep="	", header=None))
    
            # specific time-dependent impulsiveness
            specImpulsSHM = pd.DataFrame(filedata[1:, 1:],
                                            index=filedata[1:, 0],
                                            columns=filedata[0, 1:])
        
        elif file.find("Fluctuation") != -1 and file.find("ECMA") != -1:
            # ECMA-418-2:2025 Fluctuation strength
            # -----------------------------------
            filedata = np.array(pd.read_csv(file, sep="	", header=None))

            # specific ECMA-418-2:2025 Fluctuation strength
            specFluctECMA = pd.DataFrame(filedata[1:, 1:],
                                            index=filedata[1:, 0],
                                            columns=filedata[0, 1:])


            # time-dependent binaural fluctuation strength
            fluctECMASHMTDep = specFluctECMA.sum(axis=1)*bandDiff0p5

            #  fluctuation strength masked for start/end skip
            fluctECMASHMTDepMask = fluctECMASHMTDep.loc[(fluctECMASHMTDep.index.values
                                                            > start_skipT).transpose()
                                                        & (fluctECMASHMTDep.index.values
                                                            < fluctECMASHMTDep.index.values.max()
                                                            - end_skipT).transpose()]

            # overall 10% exceeded binaural fluctuation strength
            fluctECMA10Ex = fluctECMASHMTDepMask.quantile(q=0.90)

            # overall 5% exceeded fluctuation strength
            fluctECMA05Ex = fluctECMASHMTDepMask.quantile(q=0.95)

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'FluctECMA10Ex'] = fluctECMA10Ex
            dataByStim.loc[renderNames[ii], 'FluctECMA05Ex'] = fluctECMA05Ex
        
        elif file.find("Fluctuation") != -1:
            # Old Sottek Hearing Model Fluctuation strength
            # -----------------------------------
            filedata = np.array(pd.read_csv(file, sep="	", header=None))

            # specific fluctuation strength
            specFluctOldSHM = pd.DataFrame(filedata[1:, 1:],
                                           index=filedata[1:, 0],
                                           columns=filedata[0, 1:])

            # time-dependent fluctuation strength
            FluctOldSHMTDep = specFluctOldSHM.sum(axis=1)*bandDiff0p5

            # mask for start/end skip
            FluctOldSHMTDepMask = FluctOldSHMTDep.loc[(FluctOldSHMTDep.index.values
                                                        > start_skipT)
                                                        & (FluctOldSHMTDep.index.values
                                                            < FluctOldSHMTDep.index.values.max()
                                                            - end_skipT)]

            # overall (90th percentile = 10% exceeded) fluctuation strength
            FluctOldSHM10Ex = FluctOldSHMTDepMask.quantile(q=0.90)

            # overall (95th percentile = 5% exceeded) fluctuation strength
            FluctOldSHM05Ex = FluctOldSHMTDepMask.quantile(q=0.95)

            # overall time-averaged fluctuation strength
            FluctOldSHMAvg = FluctOldSHMTDepMask.mean()

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'FluctOldSHM10Ex'] = FluctOldSHM10Ex
            dataByStim.loc[renderNames[ii], 'FluctOldSHM05Ex'] = FluctOldSHM05Ex
            dataByStim.loc[renderNames[ii], 'FluctOldSHMAvg'] = FluctOldSHMAvg

        # calculation section for SQM differences
        # NOTE: THIS SECTION RELIES ON THE ALPHABETIC ORDER OF THE STIMULI FILES AS
        # ORIGINALLY NAMED: EACH AMBIENT SOUND FILE PRECEDING THE CORRESPONDING
        # COMBINED STIMULI FILES. IF THE ORDERING OR FILE NAMING IS CHANGED, THIS
        # CALCULATION WILL BE INVALIDATED AND *MAY ALSO* CAUSE AN ERROR.
        # a rolling 50 ms window is applied to average the SQM values over time -
        # this is to reduce uncertainty due to imperfect time-alignment between the
        # ambient vs combined stimuli files (all are from recordings, so there will
        # be some slippage due to imperfect editing)
        if renderNames[ii].__contains__("_Baseline"):
            # NOTE: we could dropna() the first <windowT values, but these will be
            # ignored anyway in the statistical analysis, assuming start_skipT >
            # windowT

            # calculate moving average values for ambient stimulus
            if file.find("Impulsiveness") != -1:
                ambSpecImpulsSHMMovAvg = specImpulsSHM.rolling(window=int(np.ceil(sampleRateImpulsSHM*windowT))).mean()
            
            elif file.find("Fluctuation") != -1 and file.find("ECMA") != -1:
                ambSpecFluctECMAMovAvg = specFluctECMA.rolling(window=int(np.ceil((sampleRateModECMA)*windowT))).mean()
            
        elif renderNames[ii].__contains__("Street") or renderNames[ii].__contains__("Park"):
            # calculate moving average values for combined stimulus
            # impulsiveness
            if file.find("Impulsiveness") != -1:
                specImpulsSHMMovAvg = specImpulsSHM.rolling(window=int(np.ceil(sampleRateImpulsSHM*windowT))).mean()

                dImpulsSHMTDep = np.maximum(specImpulsSHMMovAvg
                                            - ambSpecImpulsSHMMovAvg, 0).sum(axis=1).to_frame(name="Monaural")
    
                # calculate aggregated difference values
                # impulsiveness masked for start/end skip
                dImpulsSHMTDepMask = dImpulsSHMTDep.loc[(dImpulsSHMTDep.index.values
                                                            > start_skipT).transpose()
                                                        & (dImpulsSHMTDep.index.values
                                                            < dImpulsSHMTDep.index.values.max()
                                                            - end_skipT).transpose()]
                
                # overall (power-averaged) impulsiveness
                dImpulsSHMPowAvg = (dImpulsSHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).max()

                # overall (mean) impulsiveness
                dImpulsSHMAvg = dImpulsSHMTDepMask.mean().iloc[0]
    
                # overall 5% exceeded impulsiveness
                dImpulsSHM05Ex = dImpulsSHMTDepMask.quantile(q=0.95).iloc[0]

                # add results to output DataFrame
                dataByStim.loc[renderNames[ii], 'dImpulsSHMPowAvg'] = dImpulsSHMPowAvg
                dataByStim.loc[renderNames[ii], 'dImpulsSHMAvg'] = dImpulsSHMAvg
                dataByStim.loc[renderNames[ii], 'dImpulsSHM05Ex'] = dImpulsSHM05Ex

            # ECMA Sottek Hearing Model fluctuation strength
            elif file.find("Fluctuation") != -1 and file.find("ECMA") != -1:
                specFluctECMAMovAvg = specFluctECMA.rolling(window=int(np.ceil((sampleRateModECMA)*windowT))).mean()
            
                # calculate differences and make negative values 0
                # fluctuation strength
                dSpecFluctECMA = np.maximum(specFluctECMAMovAvg
                                            - ambSpecFluctECMAMovAvg, 0)

                # time-dependent fluctuation strength
                dFluctECMASHMTDep = dSpecFluctECMA.sum(axis=1)*bandDiff0p5

                # fluctuation strength masked for start/end skip
                dFluctECMASHMTDepMask = dFluctECMASHMTDep.loc[(dFluctECMASHMTDep.index.values
                                                                > start_skipT).transpose()
                                                                & (dFluctECMASHMTDep.index.values
                                                                    < dFluctECMASHMTDep.index.values.max()
                                                                    - end_skipT).transpose()]

                # overall 10% exceeded fluctuation strength
                dFluctECMA10Ex = dFluctECMASHMTDepMask.quantile(q=0.90)

                # overall 5% exceeded fluctuation strength
                dFluctECMA05Ex = dFluctECMASHMTDepMask.quantile(q=0.95)

                # add results to output DataFrame
                dataByStim.loc[renderNames[ii], 'dFluctECMA10Ex'] = dFluctECMA10Ex
                dataByStim.loc[renderNames[ii], 'dFluctECMA05Ex'] = dFluctECMA05Ex
            
    # if-else branch for non-specific sound quality
    else:
        # calculation section for Sottek Hearing Model impulsiveness
        if file.find("Impulsiveness") != -1:
        
            # read in time-dependent impulsiveness
            impulsSHMTDep = pd.read_csv(file, sep="	", header=None, index_col=0)
            
            # mask for start/end skip and 0 values
            impulsSHMTDepMask = impulsSHMTDep.loc[(impulsSHMTDep.index.values
                                                   > start_skipT).transpose()
                                                  & (impulsSHMTDep.index.values
                                                     < impulsSHMTDep.index.values.max()
                                                     - end_skipT).transpose()]

            # overall (power-averaged) impulsiveness
            impulsSHMPowAvg = (impulsSHMTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall (mean) impulsiveness
            impulsSHMAvg = impulsSHMTDepMask.mean().iloc[0]

            # overall (5% exceeded) impulsiveness
            impulsSHM05Ex = impulsSHMTDepMask.quantile(q=0.95).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'ImpulsSHMPowAvg'] = impulsSHMPowAvg
            dataByStim.loc[renderNames[ii], 'ImpulsSHMAvg'] = impulsSHMAvg
            dataByStim.loc[renderNames[ii], 'ImpulsSHM05Ex'] = impulsSHM05Ex

        elif file.find("Loudness") != -1 and file.find("ISO 532-1") != -1:

            # Calculate overall ISO 532-1 loudness from time-varing loudness
            loudISO1TDep = pd.read_csv(file, sep="	", header=None, index_col=0)

            # mask for start/end skip and 0 values
            loudISO1TDepMask = loudISO1TDep.loc[(loudISO1TDep.index.values
                                                 > start_skipT).transpose()
                                                & (loudISO1TDep.index.values
                                                   < loudISO1TDep.index.values.max()
                                                   - end_skipT).transpose()]

            # overall (5% exceeded = 95th percentile) loudness
            loudISO105Ex = loudISO1TDepMask.quantile(q=0.95).iloc[0]

            # overall (power-averaged) loudness
            loudISO1PowAvg = (loudISO1TDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]
        
            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'LoudISO105Ex'] = loudISO105Ex
            dataByStim.loc[renderNames[ii], 'LoudISO1PowAvg'] = loudISO1PowAvg

        elif file.find("Loudness") != -1 and file.find("ISO 532-3") != -1:

            # Calculate overall ISO 532-3 loudness from time-varing loudness
            loudISO3TDep = pd.read_csv(file, sep="	", header=None, index_col=0)

            # mask for start/end skip
            loudISO3TDepMask = loudISO3TDep.loc[(loudISO3TDep.index.values
                                                 > start_skipT)
                                                & (loudISO3TDep.index.values
                                                   < loudISO3TDep.index.values.max()
                                                    - end_skipT)]

            # overall (power-averaged) loudness
            loudISO3PowAvg = (loudISO3TDepMask.pow(1/np.log10(2)).mean(axis=0)**np.log10(2)).iloc[0]

            # overall 5% exceeded loudness
            loudISO305Ex = loudISO3TDepMask.quantile(q=0.95).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'LoudISO3PowAvg'] = loudISO3PowAvg
            dataByStim.loc[renderNames[ii], 'LoudISO305Ex'] = loudISO305Ex
        
        
        elif file.find("Sharpness") != -1 and file.find("DIN") != -1:

            # Calculate overall DIN 45692 sharpness from time-dependent
            # sharpness
            sharpDINTDep = pd.read_csv(file, sep="	", header=None, index_col=0)

            # mask for start/end skip
            sharpDINTDepMask = sharpDINTDep.loc[(sharpDINTDep.index.values
                                                    > start_skipT)
                                                    & (sharpDINTDep.index.values
                                                        < sharpDINTDep.index.values.max()
                                                        - end_skipT)]

            # overall (power-averaged) sharpness
            sharpDINPowAvg = (sharpDINTDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall 5% exceeded sharpness
            sharpDIN05Ex = sharpDINTDepMask.quantile(q=0.95).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'SharpDINPowAvg'] = sharpDINPowAvg
            dataByStim.loc[renderNames[ii], 'SharpDIN05Ex'] = sharpDIN05Ex
                
        elif file.find("Sharpness") != -1 and file.find("von Bismarck") != -1:

            # Calculate overall von Bismarck | ISO 532-1 sharpness from
            # time-dependent sharpness
            sharpvBISO1TDep = pd.read_csv(file, sep="	", header=None, index_col=0)

            # mask for start/end skip
            sharpvBISO1TDepMask = sharpvBISO1TDep.loc[(sharpvBISO1TDep.index.values
                                                       > start_skipT)
                                                      & (sharpvBISO1TDep.index.values
                                                         < sharpvBISO1TDep.index.values.max()    
                                                         - end_skipT)]

            # overall (power-averaged) sharpness
            sharpvBISO1PowAvg = (sharpvBISO1TDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall 5% exceeded sharpness
            sharpvBISO105Ex = sharpvBISO1TDepMask.quantile(q=0.95).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'SharpvBISO1PowAvg'] = sharpvBISO1PowAvg
            dataByStim.loc[renderNames[ii], 'SharpvBISO105Ex'] = sharpvBISO105Ex
        
        elif file.find("Sharpness") != -1 and file.find("ISO 532-3") != -1 and file.find("Aures") != -1:

            # Calculate overall Aures+ISO532-3 sharpness from time-dependent
            # sharpness
            sharpAISO3TDep = pd.read_csv(file, sep="	", header=None, index_col=0)

            # mask for start/end skip
            sharpAISO3TDepMask = sharpAISO3TDep.loc[(sharpAISO3TDep.index.values
                                                    > start_skipT)
                                                    & (sharpAISO3TDep.index.values
                                                    < sharpAISO3TDep.index.values.max()
                                                    - end_skipT)]

            # overall (power-averaged) sharpness
            sharpAISO3PowAvg = (sharpAISO3TDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall 5% exceeded sharpness
            sharpAISO305Ex = sharpAISO3TDepMask.quantile(q=0.95).iloc[0]

            # add results to output DataFrame    
            dataByStim.loc[renderNames[ii], 'SharpAurISO3PowAvg'] = sharpAISO3PowAvg
            dataByStim.loc[renderNames[ii], 'SharpAurISO305Ex'] = sharpAISO305Ex
        
        elif file.find("Sharpness") != -1 and file.find("ISO 532-1") != -1 and file.find("Aures") != -1:

            # Calculate overall Aures+ISO532-1 sharpness from time-dependent
            # sharpness
            sharpAISO1TDep = pd.read_csv(file, sep="	", header=None, index_col=0)

            # mask for start/end skip
            sharpAISO1TDepMask = sharpAISO1TDep.loc[(sharpAISO1TDep.index.values
                                                    > start_skipT)
                                                    & (sharpAISO1TDep.index.values
                                                    < sharpAISO1TDep.index.values.max()
                                                    - end_skipT)]

            # overall (power-averaged) sharpness
            sharpAISO1PowAvg = (sharpAISO1TDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall 5% exceeded sharpness
            sharpAISO105Ex = sharpAISO1TDepMask.quantile(q=0.95).iloc[0]

            # overall median sharpness
            sharpAISO1Med = sharpAISO1TDepMask.median().iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'SharpAurISO1Med'] = sharpAISO1Med
            dataByStim.loc[renderNames[ii], 'SharpAurISO1PowAvg'] = sharpAISO1PowAvg
            dataByStim.loc[renderNames[ii], 'SharpAurISO105Ex'] = sharpAISO105Ex

        # calculation section for SQM differences
        # NOTE: THIS SECTION RELIES ON THE ALPHABETIC ORDER OF THE STIMULI FILES AS
        # ORIGINALLY NAMED: EACH AMBIENT SOUND FILE PRECEDING THE CORRESPONDING
        # COMBINED STIMULI FILES. IF THE ORDERING OR FILE NAMING IS CHANGED, THIS
        # CALCULATION WILL BE INVALIDATED AND *MAY ALSO* CAUSE AN ERROR.
        # a rolling 50 ms window is applied to average the SQM values over time -
        # this is to reduce uncertainty due to imperfect time-alignment between the
        # ambient vs combined stimuli files (all are from recordings, so there will
        # be some slippage due to imperfect editing)
        if renderNames[ii].__contains__("_Baseline") and (file.find("Sharpness") != -1
                                                          and file.find("ISO 532-3") != -1
                                                          and file.find("Aures") != -1):
            # NOTE: we could dropna() the first <windowT values, but these will be
            # ignored anyway in the statistical analysis, assuming start_skipT >
            # windowT

            # calculate moving average values for ambient stimulus
            ambSharpAISO3TDepMovAvg = sharpAISO3TDep.rolling(window=int(np.ceil(sampleRateLoudISO3*windowT))).mean()

        elif (renderNames[ii].__contains__("Street")
              or renderNames[ii].__contains__("Park")) and (file.find("Sharpness") != -1
                                                            and file.find("ISO 532-3") != -1
                                                            and file.find("Aures") != -1):
            # calculate moving average values for combined stimulus
            sharpAISO3TDepMovAvg = sharpAISO3TDep.rolling(window=int(np.ceil(sampleRateLoudISO3*windowT))).mean()

            # calculate differences and make negative values 0
            dSharpAISO3TDep = np.maximum(sharpAISO3TDepMovAvg
                                         - ambSharpAISO3TDepMovAvg, 0)

            # calculate aggregated difference values
            # sharpness masked for start/end skip
            dSharpAISO3TDepMask = dSharpAISO3TDep.loc[(dSharpAISO3TDep.index.values
                                                       > start_skipT).transpose()
                                                      & (dSharpAISO3TDep.index.values
                                                         < dSharpAISO3TDep.index.values.max()
                                                         - end_skipT).transpose()]

            # overall (power-averaged) sharpness
            dSharpAISO3PowAvg = (dSharpAISO3TDepMask.pow(1/np.log10(2)).mean()**np.log10(2)).iloc[0]

            # overall 5% exceeded sharpness
            dSharpAISO305Ex = dSharpAISO3TDepMask.quantile(q=0.95).iloc[0]

            # add results to output DataFrame
            dataByStim.loc[renderNames[ii], 'dSharpAurISO3PowAvg'] = dSharpAISO3PowAvg
            dataByStim.loc[renderNames[ii], 'dSharpAurISO305Ex'] = dSharpAISO305Ex

# end of for loop over ArtemiS SQM files

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Psychoacoustic annoyance metrics single values
# ----------------------------------------------

dataByStim['PsychAnnoyWidmann'] = psych_annoy.widmannPA(dataByStim.loc[:, 'LoudISO105Ex'],
                                                        dataByStim.loc[:, 'SharpDIN05Ex'],
                                                        dataByStim.loc[:, 'RoughFZ05Ex'],
                                                        dataByStim.loc[:, 'FluctFZ05Ex'])

dataByStim['PsychAnnoyMore'] = psych_annoy.morePA(dataByStim.loc[:, 'LoudISO105Ex'],
                                                  dataByStim.loc[:, 'SharpvBISO105Ex'],
                                                  dataByStim.loc[:, 'RoughFZ05Ex'],
                                                  dataByStim.loc[:, 'FluctFZ05Ex'],
                                                  dataByStim.loc[:, 'TonalAur05Ex'])

dataByStim['PsychAnnoyDi'] = psych_annoy.diPA(dataByStim.loc[:, 'LoudISO105Ex'],
                                              dataByStim.loc[:, 'SharpvBISO105Ex'],
                                              dataByStim.loc[:, 'RoughECMAAvg'],
                                              dataByStim.loc[:, 'FluctOldSHMAvg'],
                                              dataByStim.loc[:, 'TonalAurAvg'])

dataByStim['PsychAnnoyTorija'] = psych_annoy.torijaPA(dataByStim.loc[:, 'LoudISO105Ex'],
                                                      dataByStim.loc[:, 'SharpDIN05Ex'],
                                                      dataByStim.loc[:, 'RoughECMA05Ex'],
                                                      dataByStim.loc[:, 'FluctECMA05Ex'],
                                                      dataByStim.loc[:, 'TonalECMA05Ex'],
                                                      dataByStim.loc[:, 'ImpulsSHM05Ex'])

dataByStim['PsychAnnoyWillemsen'] = psych_annoy.willemsenPA(dataByStim.loc[:, 'LoudISO105Ex'],
                                                            dataByStim.loc[:, 'SharpAurISO1Med'],
                                                            dataByStim.loc[:, 'RoughECMA05Ex'],
                                                            dataByStim.loc[:, 'ImpulsLoudWZAvg'])

dataByStim['PsychAnnoyBoucher'] = psych_annoy.boucherPA(dataByStim.loc[:, 'LoudISO105Ex'],
                                                        dataByStim.loc[:, 'SharpDIN05Ex'],
                                                        dataByStim.loc[:, 'RoughECMAAvg'],
                                                        dataByStim.loc[:, 'FluctOldSHMAvg'],
                                                        dataByStim.loc[:, 'TonalECMAAvg'])


# for rows in dataByStim.index that do not include "Baseline", fill in the dPsychAnnoy... columns by
# subtracting the row values that include with the same ambient environment in the AmbientEnv field
# from the corresponding "Baseline" row values 
for row in dataByStim.index:
    if not row.__contains__("Baseline") and not row.__contains__("Background"):
        if dataByStim.loc[row, 'AmbientRef'] == "UAS only":
            continue
        else:
            ambient_ref = dataByStim.loc[row, 'AmbientRef']
            baseline_row = dataByStim.loc[(dataByStim.AmbientRef == ambient_ref)
                                            & (dataByStim.index.str.contains("Baseline")), :]
            for col in ['PsychAnnoyWidmann',
                        'PsychAnnoyMore',
                        'PsychAnnoyDi',
                        'PsychAnnoyTorija',
                        'PsychAnnoyWillemsen',
                        'PsychAnnoyBoucher']:
                dataByStim.loc[row, 'd' + col] = (dataByStim.loc[row, col]
                                                  - baseline_row[col].iloc[0])

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# add UAS-only and ambient-only data to combined stimuli
# ------------------------------------------------------
indicesDiffPsycho = ["PartLoudSHMPowAvg",
                     "PartTonLdSHMPowAvg",
                     "PartSharpAurSHMPowAvg",
                     "PartSharpAurSHM05Ex",
                     "PartSharpWidSHMPowAvg",
                     "PartSharpWidSHM05Ex",
                     "PartSharpvBSHMPowAvg",
                     "PartSharpvBSHM05Ex",
                     "PartTonShpAurSHMPowAvg",
                     "PartTonShpAurSHM05Ex",
                     "PartTonShpWidSHMPowAvg",
                     "PartTonShpWidSHM05Ex",
                     "PartTonShpvBSHMPowAvg",
                     "PartTonShpvBSHM05Ex",
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

# list ambient references excluding "UAS only" and NaNs
ambientRefs = dataByStim.AmbientRef.unique().tolist()
ambientRefs.remove("UAS only")
ambientRefs = [ref for ref in ambientRefs if isinstance(ref, str)]

uasOnly = dataByStim.loc[dataByStim.AmbientRef == "UAS only", indices].copy()
ambOnly = dataByStim.loc[dataByStim.index.str.find("Baseline") != -1, indices].copy()
uasOnly = uasOnly.add_prefix("UAS", axis=1)
ambOnly = ambOnly.add_prefix("Amb", axis=1)

uasAmb = pd.DataFrame()
for ambRef in ambientRefs:
    # select UAS only rows for this ambient reference (except Baseline)
    maskStim = (dataByStim.AmbientRef == ambRef) & (dataByStim.index.str.find("Baseline") == -1)
    # create a logical mask of the uasOnly where the maskStim index is found in
    # the masked dataByStim index after the ambRef portion has been omitted
    maskUAS = uasOnly.index.isin(dataByStim.index[maskStim].str.replace(ambRef + "_", ""))
    uasOnlyRef = uasOnly.loc[maskUAS, :].copy()
    ambOnlyRef = ambOnly.loc[ambOnly.index.str.find(ambRef + "_") != -1, :].copy()
    # duplicate the ambOnlyRef rows to match the number of uasOnlyRef rows and add
    # the uasOnlyRef indices to the end of the ambOnlyRef index
    ambOnlyRef = ambOnlyRef.loc[ambOnlyRef.index.repeat(len(uasOnlyRef) + 1)]
    ambOnlyRef['newindex'] = ambOnlyRef.index
    ambOnlyRef.iloc[1:, -1] = [old_index.replace("_Baseline", "") + "_" + new_index for old_index, new_index in zip(ambOnlyRef.iloc[1:, :].index, uasOnlyRef.index)]
    ambOnlyRef.set_index('newindex', inplace=True)

    # rename the uasOnlyRef index to include the ambRef portion (this must happen after the ambOnlyRef reindexing)
    uasOnlyRef['newindex'] = [ambRef + "_" + file for file in list(uasOnlyRef.index)]
    uasOnlyRef.set_index('newindex', inplace=True)
    
    # add a row to uasOnlyRef for the Baseline condition with all NaNs
    baselineRow = pd.DataFrame(np.nan, index=[ambRef + "_Baseline"],
                               columns=uasOnlyRef.columns)
    uasOnlyRef = pd.concat([baselineRow, uasOnlyRef], axis=0)

    # combine the uasOnlyRef and ambOnlyRef data
    uasAmbTemp = pd.concat([uasOnlyRef, ambOnlyRef], axis=1)

    # concatenate results
    uasAmb = pd.concat([uasAmb, uasAmbTemp], axis=0)


# insert zeros for No UAS stimuli UAS SQMs
maskNoUASStims = uasAmb.index.str.contains("Baseline")
indicesUASPsycho = ["UAS" + index for index in indicesPsycho if index not in indicesDiffPsycho]
uasAmb.loc[maskNoUASStims, indicesUASPsycho] = 0

# add UAS and ambient only data to dataByStim
dataByStim = dataByStim.merge(uasAmb.astype(float), how='outer',
                              left_index=True, right_index=True)

# insert zeros for SQM differences in baseline stimuli
dataByStim.loc[dataByStim.index.str.contains("Baseline"), indicesDiffPsycho] = 0

# calculate level differences and ratios between UAS and ambient
dataByStim['LAeqLAF90diff'] = dataByStim['UASLAeqMaxLR'] - dataByStim['AmbLAF90ExMaxLR']
dataByStim['LAeqLAF10diff'] = dataByStim['UASLAeqMaxLR'] - dataByStim['AmbLAF10ExMaxLR']
dataByStim['LAF10LAF10diff'] = dataByStim['UASLAF10ExMaxLR'] - dataByStim['AmbLAF10ExMaxLR']
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

indicesLevelDiffs = ['LAeqLAF90diff', 'LAeqLAF10diff', 'LAF10LAF10diff',
                     'LASmaxLAF90diff', 'LASmaxLAF50diff',
                     'LASmaxLAeqdiff', 'LAELAF90diff', 'LAELAF50diff',
                     'LAELAeqdiff', 'PNLMLAeqdiff', 'PNLTMLAeqdiff',
                     'EPNLLAeqdiff', 'PNLMLAF90diff', 'PNLTMLAF90diff',
                     'EPNLLAF90diff', 'PNLMLAF50diff', 'PNLTMLAF50diff',
                     'EPNLLAF50diff']

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Detection-discounted UAS sound levels LAeq and LAE
# --------------------------------------------------

dataByStim['UASDisc0p5LAeqMaxLR'] = dataByStim['UASLAeqMaxLR'] - dataByStim['Detect0p5dBADiscMaxLR']
dataByStim['UASDisc0p5LAEMaxLR'] = dataByStim['UASLAEMaxLR'] - dataByStim['Detect0p5dBADiscMaxLR']
dataByStim['UASDisc0p1LAeqMaxLR'] = dataByStim['UASLAeqMaxLR'] - dataByStim['Detect0p1dBADiscMaxLR']
dataByStim['UASDisc0p1LAEMaxLR'] = dataByStim['UASLAEMaxLR'] - dataByStim['Detect0p1dBADiscMaxLR']

indicesDetect = indicesDetect + ['UASDisc0p5LAeqMaxLR', 'UASDisc0p5LAEMaxLR',
                                 'UASDisc0p1LAeqMaxLR', 'UASDisc0p1LAEMaxLR']


# %%%%%%%%%%%%%%%%%%%%%%%
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

# move detectability metrics to after level difference metrics
indicesDetect.reverse()
cols = list(dataByStim.columns)
for col in indicesDetect:
    cols.insert(cols.index(indicesLevelDiffs[-1]) + 1,
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

# move partial loudness metrics to after level difference metrics
cols = list(dataByStim.columns)
cols.insert(cols.index(indicesLevelDiffs[-1]) + 1,
            cols.pop(cols.index('PartLoudSHMPowAvg')))
cols.insert(cols.index('PartLoudSHMPowAvg') + 1,
            cols.pop(cols.index('PartTonLdSHMPowAvg')))
# for each partial sharpness and partial tonal sharpness metric, move these
# to after partial loudness metrics
partSharpMetrics = ['PartSharpAurSHMPowAvg',
                    'PartSharpAurSHM05Ex',
                    'PartSharpWidSHMPowAvg',
                    'PartSharpWidSHM05Ex',
                    'PartSharpvBSHMPowAvg',
                    'PartSharpvBSHM05Ex',
                    'PartTonShpAurSHMPowAvg',
                    'PartTonShpAurSHM05Ex',
                    'PartTonShpWidSHMPowAvg',
                    'PartTonShpWidSHM05Ex',
                    'PartTonShpvBSHMPowAvg',
                    'PartTonShpvBSHM05Ex']
for ii, col in enumerate(partSharpMetrics):
    cols.insert(cols.index('PartTonLdSHMPowAvg') + ii + 1,
                cols.pop(cols.index(col)))
dataByStim = dataByStim[cols]

# set all partial metrics to 0 for Baseline (no UAS) stimuli
maskBaselineStims = dataByStim.index.str.contains("Baseline")
partialMetrics = ['PartLoudSHMPowAvg',
                  'PartTonLdSHMPowAvg'] + partSharpMetrics
dataByStim.loc[maskBaselineStims, partialMetrics] = 0

# %%%%%%%%%%%%%
# Response data
# -------------

# open csv file selection dialog and assign filepath
# PROJECT NOTE: the results files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%202/Test_files/Response_data/Compiled?csf=1&web=1&e=pl8nrz
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

fileExts = "*.csv"
filepath = QFileDialog.getOpenFileName(filter=fileExts,
                                       caption=r"Select test end response data file in '03 Experiment\Experiment 2\Test_files\Response_data\Compiled'")[0]

# read in data and add column indicating the stimulus recording file
testResponses = pd.read_csv(filepath, header=0)

# sort stimulus names
stimSorted = np.sort(testResponses['stimulus'].unique())

# initialise DataFrames for loop over stimuli
testData = pd.DataFrame(index=stimSorted)

# loop over stimuli recording names, extract corresponding response data and
# tranpose data with recording as index
# calculate aggregate statistics
for ii, file in enumerate(stimSorted):
    if ii == 0:
        print("Processing results...\n")
    print(file)

    # loop for each response type, extract individual responses,
    # and calculate median and mean aggregations
    for response in ["Pleasantness", "dPleasantness",
                     "Eventfulness", "dEventfulness",
                     "Annoyance", "dAnnoyance"]:
        responseData = testResponses.loc[testResponses['stimulus'] == file,
                                         ['participant', response]]
        columns = [response + "_" + str(ID) for ID in responseData['participant']]
        responseData = pd.DataFrame(data=np.array(responseData[response]),
                                    index=columns, columns=[file]).transpose()
        
        # if responseData is not all NaN (no responses), calculate median and mean
        if not np.all(responseData.isna()):
            responseAgg = pd.DataFrame(data=np.vstack([np.nanpercentile(responseData.values,
                                                                        q=50, axis=1,
                                                                        method='median_unbiased')[0],
                                                    np.nanmean(responseData.values, axis=1)[0]]),
                                       index=[response + 'Median',
                                              response + 'Mean'],
                                       columns=[file]).transpose()
        else:
            responseAgg = pd.DataFrame(data=np.vstack([np.nan, np.nan]),
                                       index=[response + 'Median',
                                              response + 'Mean'],
                                       columns=[file]).transpose()

        # add to testData DataFrame
        if ii == 0:
            testData = testData.join(responseAgg, how='outer')
            testData = testData.join(responseData, how='outer')
        else:
            testData.loc[file, response + 'Median'] = responseAgg.loc[file, response + 'Median']
            testData.loc[file, response + 'Mean'] = responseAgg.loc[file, response + 'Mean']
            
            for col in columns:
                testData.loc[file, col] = responseData.loc[file, col]
        

    for response in ["HighlyAnnoyed", "dHighlyAnnoyed"]:
        responseData = testResponses.loc[testResponses['stimulus'] == file,
                                         ['participant', response]]
        columns = [response + "_" + str(ID) for ID in responseData['participant']]
        responseData = pd.DataFrame(data=np.array(responseData[response]),
                                    index=columns, columns=[file]).transpose()
        
        # if responseData is not all NaN (no responses), calculate total and proportion
        if not np.all(responseData.isna()):
            responseAgg = pd.DataFrame(data=np.vstack([np.nansum(responseData.values,
                                                                 axis=1)[0],
                                                       np.nanmean(responseData.values,
                                                                  axis=1)[0]]),
                                       index=[response + 'Total',
                                              response + 'Prop'],
                                       columns=[file]).transpose()
        else:
            responseAgg = pd.DataFrame(data=np.vstack([np.nan, np.nan]),
                                       index=[response + 'Total',
                                              response + 'Prop'],
                                       columns=[file]).transpose()

        # add to testData DataFrame
        if ii == 0:
            testData = testData.join(responseAgg, how='outer')
            testData = testData.join(responseData, how='outer')
        else:
            testData.loc[file, response + 'Total'] = responseAgg.loc[file, response + 'Total']
            testData.loc[file, response + 'Prop'] = responseAgg.loc[file, response + 'Prop']
            
            for col in columns:
                testData.loc[file, col] = responseData.loc[file, col]

# move all the individual response columns to the front of the DataFrame
# and group by response type
indivCols = [col for col in testData.columns
             if not (col.endswith('Median')
                     or col.endswith('Mean')
                     or col.endswith('Total')
                     or col.endswith('Prop'))]
indivCols.sort()
cols = list(testData.columns)
cols.sort()
for col in reversed(indivCols):
    cols.insert(0, cols.pop(cols.index(col)))
testData = testData[cols]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# combining left and right data
# -----------------------------

# form a dataByStimCombi DataFrame that uses the artihmetic mean of 
# every parameter across the left or right versions of each stimulus
# first make a dataByStimL by filtering dataByStim for UASStart = Left
dataByStimL = dataByStim.loc[dataByStim.UASStart == "Left", :].copy()
dataByStimL.index = dataByStimL.index.str.replace("_left", "")
dataByStimL.drop(columns=['HATSRecFiles', 'MA220MicRecFiles', 'UASStart'], inplace=True)
# now do right
dataByStimR = dataByStim.loc[dataByStim.UASStart == "Right", :].copy()
dataByStimR.index = dataByStimR.index.str.replace("_right", "")
dataByStimR.drop(columns=['HATSRecFiles', 'MA220MicRecFiles', 'UASStart'], inplace=True)
# now calculate the element wise mean across the two dataFrames, across columns from
# 'LAeqMaxLR' to the end 
dataByStimCombi = dataByStimL.copy()
dataByStimCombi.loc[:, 'LAeqMaxLR':] = 0.5*(dataByStimL.loc[:, 'LAeqMaxLR':].values
                                            + dataByStimR.loc[:, 'LAeqMaxLR':].values)
# assign new stimID values to the combined stimuli by starting the index from the maximum
# stimID value in dataByStim and adding 1 for each row in dataByStimCombi
maxStimID = int(dataByStim['StimID'].max())
dataByStimCombi['StimID'] = range(maxStimID + 1, maxStimID + 1 + len(dataByStimCombi))

# now merge this with the remaining rows from dataByStim that do not include "Left" or "Right" in the index
dataByStimCombi = pd.concat([dataByStim.loc[~dataByStim.index.str.contains("left|right"), :].drop(columns=['HATSRecFiles', 'MA220MicRecFiles']),
                             dataByStimCombi], axis=0)
# resort the index
dataByStimCombi.sort_index(inplace=True)
# drop the UASStart column from dataByStimCombi as this is no longer relevant
dataByStimCombi.drop(columns=['UASStart'], inplace=True)

# also combine the testData into a combined dataFrame
testDataL = testData.loc[testData.index.str.contains("left"), :].copy()
testDataL.index = testDataL.index.str.replace("_left", "")
testDataR = testData.loc[testData.index.str.contains("right"), :].copy()
testDataR.index = testDataR.index.str.replace("_right", "")
testDataCombi = testDataL.fillna(testDataR)
testDataCombi = pd.concat([testData.loc[~testData.index.str.contains("left|right"), :],
                           testDataCombi], axis=0)
testDataCombi.sort_index(inplace=True)


# %%%%%%%%%%%%%%%%%%
# Merge the datasets
# ------------------

# merge testData into dataByStim, matching on the index
dataByStim = dataByStim.merge(testData, how='outer', left_index=True,
                              right_index=True)

# merge testDataCombi into dataByStimCombi, matching on the index
dataByStimCombi = dataByStimCombi.merge(testDataCombi, how='outer', left_index=True,
                                        right_index=True)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# open questionnaire responses dataset
# ------------------------------------

# open csv file selection dialog and assign filepath
# PROJECT NOTE: the response files are stored in
# https://testlivesalfordac.sharepoint.com/:f:/r/sites/REFMAP/Shared%20Documents/General/03%20Experiment/Experiment%202/Test_files/Questionnaire?csf=1&web=1&e=W4YAdI
# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

fileExts = "*.csv"
filepath = QFileDialog.getOpenFileName(filter=fileExts,
                                       caption=r"Select test questionnaire response data file in '03 Experiment\Experiment 2\Test_files\Questionnaire'")[0]

questResponses = pd.read_csv(filepath, header=0)
# convert questResponses Exp1ID column to int64 type
questResponses['Exp1ID'] = questResponses['Exp1ID'].astype(pd.Int64Dtype())
questResponses['AAMExperience'] = questResponses['AAMExperience'].fillna("None")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Prepare outputs for saving to file
# ----------------------------------

# separate 'by stimulus' output into test data and auxiliary data by checking
# for numeric StimID values (test stimuli only have numeric StimID values)
dataByStimTest = dataByStim.loc[~np.isnan(dataByStim['StimID']), :]
dataByStimTestCombi = dataByStimCombi.loc[~np.isnan(dataByStimCombi['StimID']), :]

# select auxiliary data (non-test stimuli) by indexes with background, and the UAS only stimuli
dataByStimAux = dataByStim.loc[(dataByStim.index.str.contains("Background"))
                               | (dataByStim['AmbientRef'] == "UAS only"), :].dropna(axis=1, how='all')
dataByStimAuxCombi = dataByStimCombi.loc[(dataByStimCombi.index.str.contains("Background"))
                                          | (dataByStimCombi['AmbientRef'] == "UAS only"), :].dropna(axis=1, how='all')

# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 

outFilePath = QFileDialog.getExistingDirectory(caption="Choose output folder to save processed files in '03 Experiment\Experiment 2\Analysis\PostProcess'")

dataByStim.to_csv(os.path.join(outFilePath,
                               "refmap_listest2_alldata_ByStim.csv"))

dataByStimTest.to_csv(os.path.join(outFilePath,
                                   "refmap_listest2_testdata_ByStim.csv"))

dataByStimAux.to_csv(os.path.join(outFilePath,
                                  "refmap_listest2_auxdata.csv"))

dataByStimCombi.to_csv(os.path.join(outFilePath,
                                    "refmap_listest2_alldata_ByStimCombi.csv"))

dataByStimTestCombi.to_csv(os.path.join(outFilePath,
                                        "refmap_listest2_testdata_ByStimCombi.csv"))

dataByStimAuxCombi.to_csv(os.path.join(outFilePath,
                                       "refmap_listest2_auxdata_ByStimCombi.csv"))

# merge response and stimuli data into 'by participant' test datasets, and
# save to file

testDataBySubj = pd.merge(left=testResponses.drop(columns=['ambientRef', 'sourceType',
                                                           'sourceMode', 'sourceProximity',
                                                           'sourceStart', 'sourceEvents', 'sourceInterval']),
                          right=dataByStimTest.loc[:, :dataByStimTest.columns[dataByStimTest.columns.get_loc('Annoyance_1') - 1]],
                          how='outer', left_on='stimulus', right_index=True)

testDataBySubj = pd.merge(left=testDataBySubj,
                          right=questResponses, how='left',
                          left_on='participant', right_on='ParticipantID')

# rename participant column to ID
testDataBySubj.rename(columns={'participant': 'ID', 'trial': 'Trial', 'stimulus': 'Stimulus'}, inplace=True)
testDataBySubj.drop(columns=['ParticipantID'], inplace=True)
testDataBySubj.sort_values(by='ID', axis=0, inplace=True)

testDataBySubj.to_csv(os.path.join(outFilePath,
                                  "refmap_listest2_testdata_BySubj.csv"),
                      index=False)

# form wide format datasets for each outcome
testAnnoyDataBySubjWide = testDataBySubj.pivot(index='ID',
                                               columns='Stimulus',
                                               values='Annoyance')
testPleasantDataBySubjWide = testDataBySubj.pivot(index='ID',
                                                  columns='Stimulus',
                                                  values='Pleasantness')
testEventfulDataBySubjWide = testDataBySubj.pivot(index='ID',
                                                  columns='Stimulus',
                                                  values='Eventfulness')
testdAnnoyDataBySubjWide = testDataBySubj.pivot(index='ID',
                                                columns='Stimulus',
                                                values='dAnnoyance')
testdPleasantDataBySubjWide = testDataBySubj.pivot(index='ID',
                                                   columns='Stimulus',
                                                   values='dPleasantness')
testdEventfulDataBySubjWide = testDataBySubj.pivot(index='ID',
                                                   columns='Stimulus',
                                                   values='dEventfulness')

# merge with participant info using ID and ParticipantID, dropping ParticipantID
testAnnoyDataBySubjWide = testAnnoyDataBySubjWide.merge(questResponses,
                                                        left_on='ID', right_on='ParticipantID',
                                                        how='left').rename(columns={'ParticipantID': 'ID'})
# move ID column to front
cols = list(testAnnoyDataBySubjWide.columns)
cols.insert(0, cols.pop(cols.index('ID')))
testAnnoyDataBySubjWide = testAnnoyDataBySubjWide[cols]

testPleasantDataBySubjWide = testPleasantDataBySubjWide.merge(questResponses,
                                                              left_on='ID', right_on='ParticipantID',
                                                              how='left').rename(columns={'ParticipantID': 'ID'})
cols = list(testPleasantDataBySubjWide.columns)
cols.insert(0, cols.pop(cols.index('ID')))
testPleasantDataBySubjWide = testPleasantDataBySubjWide[cols]

testEventfulDataBySubjWide = testEventfulDataBySubjWide.merge(questResponses,
                                                              left_on='ID', right_on='ParticipantID',
                                                              how='left').rename(columns={'ParticipantID': 'ID'})
cols = list(testEventfulDataBySubjWide.columns)
cols.insert(0, cols.pop(cols.index('ID')))
testEventfulDataBySubjWide = testEventfulDataBySubjWide[cols]

testdAnnoyDataBySubjWide = testdAnnoyDataBySubjWide.merge(questResponses,
                                                          left_on='ID', right_on='ParticipantID',
                                                          how='left').rename(columns={'ParticipantID': 'ID'})
cols = list(testdAnnoyDataBySubjWide.columns)
cols.insert(0, cols.pop(cols.index('ID')))
testdAnnoyDataBySubjWide = testdAnnoyDataBySubjWide[cols]

testdPleasantDataBySubjWide = testdPleasantDataBySubjWide.merge(questResponses,
                                                                left_on='ID', right_on='ParticipantID',
                                                                how='left').rename(columns={'ParticipantID': 'ID'})
cols = list(testdPleasantDataBySubjWide.columns)
cols.insert(0, cols.pop(cols.index('ID')))
testdPleasantDataBySubjWide = testdPleasantDataBySubjWide[cols]

testdEventfulDataBySubjWide = testdEventfulDataBySubjWide.merge(questResponses,
                                                                left_on='ID', right_on='ParticipantID',
                                                                how='left').rename(columns={'ParticipantID': 'ID'})
cols = list(testdEventfulDataBySubjWide.columns)
cols.insert(0, cols.pop(cols.index('ID')))
testdEventfulDataBySubjWide = testdEventfulDataBySubjWide[cols]


# save wide format datasets to file
testAnnoyDataBySubjWide.to_csv(os.path.join(outFilePath,
                                            "refmap_listest2_AnnoyDataBySubjWide.csv"),
                               index=False)

testPleasantDataBySubjWide.to_csv(os.path.join(outFilePath,
                                               "refmap_listest2_PleasantDataBySubjWide.csv"),
                                  index=False)

testEventfulDataBySubjWide.to_csv(os.path.join(outFilePath,
                                               "refmap_listest2_EventfulDataBySubjWide.csv"),
                                  index=False)

testdAnnoyDataBySubjWide.to_csv(os.path.join(outFilePath,
                                              "refmap_listest2_dAnnoyDataBySubjWide.csv"),
                                index=False)

testdPleasantDataBySubjWide.to_csv(os.path.join(outFilePath,
                                                "refmap_listest2_dPleasantDataBySubjWide.csv"),
                                   index=False)

testdEventfulDataBySubjWide.to_csv(os.path.join(outFilePath,
                                               "refmap_listest2_dEventfulDataBySubjWide.csv"),
                                   index=False)

# repeat for combined dataset
testDataCombiBySubj = pd.merge(left=testResponses.drop(columns=['ambientRef', 'sourceType',
                                                                'sourceMode', 'sourceProximity',
                                                                'sourceStart', 'sourceEvents', 'sourceInterval']),
                               right=dataByStimCombi.loc[:, :dataByStimCombi.columns[dataByStimCombi.columns.get_loc('Annoyance_1') - 1]],
                               how='outer', left_on='stimulus', right_index=True)

testDataCombiBySubj = pd.merge(left=testDataCombiBySubj,
                               right=questResponses, how='left',
                               left_on='participant', right_on='ParticipantID')

testDataCombiBySubj.rename(columns={'participant': 'ID', 'trial': 'Trial', 'stimulus': 'Stimulus'}, inplace=True)
testDataCombiBySubj.drop(columns=['ParticipantID'], inplace=True)
testDataCombiBySubj.sort_values(by='ID', axis=0, inplace=True)

testDataCombiBySubj.to_csv(os.path.join(outFilePath,
                                        "refmap_listest2_testdataCombi_BySubj.csv"),
                           index=False)

testAnnoyDataCombiBySubjWide = testDataCombiBySubj.pivot(index='ID',
                                                         columns='Stimulus',
                                                         values='Annoyance')

testPleasantDataCombiBySubjWide = testDataCombiBySubj.pivot(index='ID',
                                                            columns='Stimulus',
                                                            values='Pleasantness')

testEventfulDataCombiBySubjWide = testDataCombiBySubj.pivot(index='ID',
                                                            columns='Stimulus',
                                                            values='Eventfulness')

testdAnnoyDataCombiBySubjWide = testDataCombiBySubj.pivot(index='ID',
                                                          columns='Stimulus',
                                                          values='dAnnoyance')

testdPleasantDataCombiBySubjWide = testDataCombiBySubj.pivot(index='ID',
                                                             columns='Stimulus',
                                                             values='dPleasantness')

testdEventfulDataCombiBySubjWide = testDataCombiBySubj.pivot(index='ID',
                                                             columns='Stimulus',
                                                             values='dEventfulness')

testAnnoyDataCombiBySubjWide = testAnnoyDataCombiBySubjWide.merge(questResponses,
                                                                  left_on='ID', right_on='ParticipantID',
                                                                  how='left').rename(columns={'ParticipantID': 'ID'})

cols = list(testAnnoyDataCombiBySubjWide.columns)
cols.insert(0, cols.pop(cols.index('ID')))
testAnnoyDataCombiBySubjWide = testAnnoyDataCombiBySubjWide[cols]

testPleasantDataCombiBySubjWide = testPleasantDataCombiBySubjWide.merge(questResponses,
                                                                        left_on='ID', right_on='ParticipantID',
                                                                        how='left').rename(columns={'ParticipantID': 'ID'})

cols = list(testPleasantDataCombiBySubjWide.columns)
cols.insert(0, cols.pop(cols.index('ID')))
testPleasantDataCombiBySubjWide = testPleasantDataCombiBySubjWide[cols]

testEventfulDataCombiBySubjWide = testEventfulDataCombiBySubjWide.merge(questResponses,
                                                                        left_on='ID', right_on='ParticipantID',
                                                                        how='left').rename(columns={'ParticipantID': 'ID'})
cols = list(testEventfulDataCombiBySubjWide.columns)
cols.insert(0, cols.pop(cols.index('ID')))
testEventfulDataCombiBySubjWide = testEventfulDataCombiBySubjWide[cols]

testdAnnoyDataCombiBySubjWide = testdAnnoyDataCombiBySubjWide.merge(questResponses,
                                                                    left_on='ID', right_on='ParticipantID',
                                                                    how='left').rename(columns={'ParticipantID': 'ID'})
cols = list(testdAnnoyDataCombiBySubjWide.columns)
cols.insert(0, cols.pop(cols.index('ID')))
testdAnnoyDataCombiBySubjWide = testdAnnoyDataCombiBySubjWide[cols]

testdPleasantDataCombiBySubjWide = testdPleasantDataCombiBySubjWide.merge(questResponses,
                                                                          left_on='ID', right_on='ParticipantID',
                                                                          how='left').rename(columns={'ParticipantID': 'ID'})
cols = list(testdPleasantDataCombiBySubjWide.columns)
cols.insert(0, cols.pop(cols.index('ID')))
testdPleasantDataCombiBySubjWide = testdPleasantDataCombiBySubjWide[cols]

testdEventfulDataCombiBySubjWide = testdEventfulDataCombiBySubjWide.merge(questResponses,
                                                                          left_on='ID', right_on='ParticipantID',
                                                                          how='left').rename(columns={'ParticipantID': 'ID'})
cols = list(testdEventfulDataCombiBySubjWide.columns)
cols.insert(0, cols.pop(cols.index('ID')))
testdEventfulDataCombiBySubjWide = testdEventfulDataCombiBySubjWide[cols]

testAnnoyDataCombiBySubjWide.to_csv(os.path.join(outFilePath,
                                                 "refmap_listest2_AnnoyDataCombiBySubjWide.csv"),
                                                 index=False)
testPleasantDataCombiBySubjWide.to_csv(os.path.join(outFilePath,
                                                    "refmap_listest2_PleasantDataCombiBySubjWide.csv"),
                                                    index=False)
testEventfulDataCombiBySubjWide.to_csv(os.path.join(outFilePath,
                                                    "refmap_listest2_EventfulDataCombiBySubjWide.csv"),
                                                    index=False)
testdAnnoyDataCombiBySubjWide.to_csv(os.path.join(outFilePath,
                                                  "refmap_listest2_dAnnoyDataCombiBySubjWide.csv"),
                                                  index=False)
testdPleasantDataCombiBySubjWide.to_csv(os.path.join(outFilePath,
                                                     "refmap_listest2_dPleasantDataCombiBySubjWide.csv"),
                                                     index=False)
testdEventfulDataCombiBySubjWide.to_csv(os.path.join(outFilePath,
                                                     "refmap_listest2_dEventfulDataCombiBySubjWide.csv"),
                                                     index=False)
