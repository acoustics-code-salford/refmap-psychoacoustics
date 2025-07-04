# -*- coding: utf-8 -*-
# %% Preamble
"""
acousticSHMLoudness.py
----------------------

Returns loudness values according to ECMA-418-2:2024 (using the Sottek Hearing
Model) for an input calibrated single mono or single stereo audio (sound
pressure) time-series signal, p.

Requirements
------------
numpy
scipy
matplotlib
tqdm
bottleneck
acoustic-toolbox
refmap-psychoacoustics (metrics.ecma418_2, dsp.filterFuncs and
                        utils.formatFuncs)

Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford

Date created: 29/05/2023
Date last modified: 02/06/2025
Python version: 3.11

Copyright statement: This file and code is part of work undertaken within
the RefMap project (www.refmap.eu), and is subject to licence as detailed
in the code repository
(https://github.com/acoustics-code-salford/refmap-psychoacoustics)

As per the licensing information, please be aware that this code is WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

"""

# %% Import block
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.fft import (fft, ifft)
from src.py.metrics.ecma418_2.acousticSHMSubs import (shmResample, shmPreProc,
                                                      shmOutMidEarFilter,
                                                      shmAuditoryFiltBank,
                                                      shmSignalSegment,
                                                      shmBasisLoudness,
                                                      shmNoiseRedLowPass)
from src.py.metrics.ecma418_2.acousticSHMTonality import acousticSHMTonality
from tqdm import tqdm
import bottleneck as bn
from src.py.dsp.filterFuncs import A_weight_T
from src.py.utils.formatFuncs import roundTrad
from acoustic_toolbox.signal import rms

# %% Module settings
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['mathtext.fontset'] = 'stixsans'

plt.rc('font', size=16)  # controls default text sizes
plt.rc('axes', titlesize=18,
       labelsize=18)  # fontsize of the axes title and x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the tick labels
plt.rc('legend', fontsize=16)  # legend fontsize
plt.rc('figure', titlesize=20)  # fontsize of the figure title


# %% acousticSHMLoudness
def acousticSHMLoudness(p, sampleRateIn, axisN=0, soundField='freeFrontal',
                        waitBar=True, outPlot=False, binaural=True):
    """
    Inputs
    ------
    p : 1D or 2D array
        the input signal as single mono or stereo audio (sound
        pressure) signals

    sampleRateIn : integer
                   the sample rate (frequency) of the input signal(s)

    axisN : integer (0 or 1, default: 0)
            the time axis along which to calculate the loudness

    soundField : keyword string (default: 'freeFrontal')
                 determines whether the 'freeFrontal' or 'diffuse' field stages
                 are applied in the outer-middle ear filter, or 'noOuter' uses
                 only the middle ear stage, or 'noEar' omits ear filtering.
                 Note: these last two options are beyond the scope of the
                 standard, but may be useful if recordings made using
                 artificial outer/middle ear are to be processed using the
                 specific recorded responses.

    waitBar : keyword string (default: True)
              determines whether a progress bar displays during processing
              (set waitBar to false for doing multi-file parallel calculations)

    outPlot : Boolean (default: False)
              flag indicating whether to generate a figure from the output
              (set outplot to false for doing multi-file parallel calculations)

    binaural : Boolean (default: True)
               flag indicating whether to output combined binaural loudness for
               stereo input signal

    Returns
    -------
    loudnessSHM : dict
                  contains the output

    loudnessSHM contains the following outputs:

    specLoudness : 2D or 3D array
                   time-dependent specific loudness for each critical band
                   arranged as [time, bands(, channels)]

    specLoudnessPowAvg : 1D or 2D array
                         time-averaged specific loudness for each critical band
                         arranged as [bands(, channels)]

    specTonalLoudness : 2D or 3D array
                        time-dependent specific tonal loudness for each
                        critical band
                        arranged as [time, bands(, channels)]
    
    specNoiseLoudness : 2D or 3D array
                        time-dependent specific noise loudness for each
                        critical band
                        arranged as [time, bands(, channels)]

    loudnessTDep : 1D or 2D array
                   time-dependent overall loudness
                   arranged as [time(, channels)]

    loudnessPowAvg : 1D or 2D array
                     time-averaged overall loudness
                     arranged as [loudness(, channels)]

    bandCentreFreqs : 1D array
                      centre frequencies corresponding with each critical band
                      rate

    timeOut : 1D array
              time (seconds) corresponding with time-dependent outputs

    soundField : string
                 identifies the soundfield type applied (the input argument
                 soundField)

    If outPlot=True, a set of plots is returned illustrating the energy
    time-averaged A-weighted sound level, the time-dependent specific and
    overall loudness, with the latter also indicating the time-aggregated
    value. A set of plots is returned for each input channel.
    
    If binaural=true, a corresponding set of outputs for the binaural
    loudness are also contained in loudnessSHM.

    Assumptions
    -----------
    The input signal is calibrated to units of acoustic pressure in Pascals
    (Pa).

    Checked by:
    Date last checked:

    """
    # %% Input checks
    # Orient input matrix
    if axisN in [0, 1]:
        if axisN == 1:
            p = p.T
    else:
        raise ValueError("Input axisN must be an integer 0 or 1")

    # Check the length of the input data (must be longer than 300 ms)
    if p.shape[0] <= 300/1000*sampleRateIn:
        raise ValueError('Input signal is too short along the specified axis to calculate loudness (must be longer than 300 ms)')

    # Check the channel number of the input data
    if p.shape[1] > 2:
        raise ValueError('Input signal comprises more than two channels')

    chansIn = p.shape[1]
    if chansIn > 1:
        chans = ["Stereo left",
                 "Stereo right"]
    else:
        chans = ["Mono"]
    # end of if branch for channel number check

    # %% Define constants

    # Signal sample rate prescribed to be 48kHz (to be used for resampling), Section 5.1.1 ECMA-418-2:2024 [r_s]
    sampleRate48k = 48e3
    # defined in Section 5.1.4.1 ECMA-418-2:2024 [deltaf(f=0)]
    deltaFreq0 = 81.9289
    # Half-overlapping Bark band centre-frequency denominator constant defined in Section 5.1.4.1 ECMA-418-2:2024
    c = 0.1618

    dz = 0.5  # critical band overlap
    # half-overlapping critical band rate scale [z]
    halfBark = np.arange(0.5, 27, dz)
    # Section 5.1.4.1 Equation 9 ECMA-418-2:2024 [F(z)]
    bandCentreFreqs = (deltaFreq0/c)*np.sinh(c*halfBark)
    # Section 5.1.4.1 Equation 10 ECMA-418-2:2024 [deltaf(z)]
    dfz = np.sqrt(deltaFreq0**2 + (c*bandCentreFreqs)**2)

    # Section 8.1.1 ECMA-418-2:2024
    weight_n = 0.5331  # Equations 113 & 114 ECMA-418-2:2024 [w_n]
    # Table 12 ECMA-418-2:2024
    a = 0.2918
    b = 0.5459

    # Output sample rate based on tonality hop sizes (Section 6.2.6
    # ECMA-418-2:2024) [r_sd]
    sampleRate1875 = 48e3/256

    # %% Signal processing

    # Input pre-processing
    # --------------------
    if sampleRateIn != sampleRate48k:  # Resample signal
        p_re, _ = shmResample(p, sampleRateIn)
    else:  # don't resample
        p_re = p
    # end of if branch for resampling

    # Calculate specific loudnesses for tonal and noise components
    # ------------------------------------------------------------

    # Obtain tonal and noise component specific loudnesses from Sections 5 & 6 ECMA-418-2:2024
    tonalitySHM = acousticSHMTonality(p_re, sampleRate48k, axisN=0,
                                      soundField=soundField,
                                      waitBar=waitBar, outPlot=False)

    specTonalLoudness = tonalitySHM['specTonalLoudness']  # [N'_tonal(l,z)]
    specNoiseLoudness = tonalitySHM['specNoiseLoudness']  # [N'_noise(l,z)]

    # Section 8.1.1 ECMA-418-2:2024
    # Weight and combine component specific loudnesses
    if waitBar:
        chanIter = tqdm(range(chansIn), desc="Channels")
    else:
        chanIter = range(chansIn)

    # pre-allocate array
    specLoudness = np.zeros(specTonalLoudness.shape)
    for chan in chanIter:
        # Equation 114 ECMA-418-2:2024 [e(z)]
        maxLoudnessFuncel = a/(np.max(specTonalLoudness[:, :, chan]
                                      + specNoiseLoudness[:, :, chan], axis=1)
                               + 1e-12) + b
        # Equation 113 ECMA-418-2:2024 [N'(l,z)]
        specLoudness[:, :, chan] = (specTonalLoudness[:, :, chan]**maxLoudnessFuncel
                                    + np.abs((weight_n*specNoiseLoudness[:,
                                                                         :,
                                                                         chan])
                                             ** maxLoudnessFuncel))**(1/maxLoudnessFuncel)
    # end of loudness for loop over channels

    if chansIn == 2 and binaural:
        # Binaural loudness
        # Section 8.1.5 ECMA-418-2:2024 Equation 118 [N'_B(l,z)]
        specLoudness[:, :, 3] = np.sqrt(np.sum(specLoudness**2, 2)/2)
        specTonalLoudness[:, :, 3] = np.sqrt(np.sum(specTonalLoudness**2, 2)/2)
        specNoiseLoudness[:, :, 3] = np.sqrt(np.sum(specNoiseLoudness**2, 2)/2)
        chansOut = 3  # set number of 'channels' to stereo plus single binaural
        chans = chans + ["Binaural"]
    else:
        chansOut = chansIn  # assign number of output channels
    # end of if branch for channels

    # Section 8.1.2 ECMA-418-2:2024
    # Time-averaged specific loudness Equation 115 [N'(z)]
    specLoudnessPowAvg = (np.sum(specLoudness[57:, :, :]**(1/np.log10(2)), 0)
                          / specLoudness[57:, :, :].shape[0])**np.log10(2)

    # Section 8.1.3 ECMA-418-2:2024
    # Time-dependent loudness Equation 116 [N(l)]
    loudnessTDep = np.sum(specLoudness**dz, 1)

    # Section 8.1.4 ECMA-418-2:2024
    # Overall loudness Equation 117 [N]
    loudnessPowAvg = (np.sum(loudnessTDep[57:, :]**(1/np.log10(2)), 0)
                      / loudnessTDep[57:, :].shape[0])**np.log10(2)

    # time (s) corresponding with results output [t]
    timeOut = np.arange(0, specLoudness.shape[0])/sampleRate1875

    # %% Output plotting

    # Plot figures
    # ------------
    if outPlot:
        if waitBar:
            plotIter = tqdm(range(chansOut), desc="Channels")
        else:
            plotIter = range(chansOut)

        for chan in plotIter:
            # Plot results
            cmap_viridis = mpl.colormaps['viridis']
            chan_lab = chans[chan]
            fig, axs = plt.subplots(nrows=2, ncols=1, figsize=[10.5, 7.5],
                                    layout='constrained')

            ax1 = axs[0]
            pmesh = ax1.pcolormesh(timeOut, bandCentreFreqs,
                                   np.swapaxes(specLoudness[:, :, chan], 0, 1),
                                   cmap=cmap_viridis,
                                   vmin=0,
                                   vmax=np.ceil(np.max(loudnessTDep[:,
                                                                    chan])*10)/10,
                                   shading='gouraud')
            ax1.set(xlim=[timeOut[1],
                          timeOut[-1] + (timeOut[1] - timeOut[0])],
                    xlabel="Time, s",
                    ylim=[bandCentreFreqs[0], bandCentreFreqs[-1]],
                    yscale='log',
                    yticks=[63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],
                    yticklabels=["63", "125", "250", "500", "1k", "2k", "4k",
                                 "8k", "16k"],
                    ylabel="Frequency, Hz")
            ax1.minorticks_off()
            cbax = ax1.inset_axes([1.05, 0, 0.05, 1])
            fig.colorbar(pmesh, ax=ax1,
                         label=(r"Specific loudness,"
                                "\n"
                                r"$\mathregular{tu_{SHM}}/\mathregular{Bark_{SHM}}$"),
                         aspect=10, cax=cbax)

            ax2 = axs[1]
            ax2.plot(timeOut, loudnessPowAvg[chan]*np.ones(timeOut.size),
                     color=cmap_viridis(33/255), linewidth=1,
                     label=("Time-" + "\n" + "average"))
            ax2.plot(timeOut, loudnessTDep[:, chan],
                     color=cmap_viridis(165/255),
                     linewidth=0.75, label=("Time-" + "\n" + "dependent"))
            ax2.set(xlim=[timeOut[0], timeOut[-1] + timeOut[1] - timeOut[0]],
                    xlabel="Time, s",
                    ylim=[0, 1.1*np.ceil(np.max(loudnessTDep[:, chan])*10)/10],
                    ylabel=(r"Loudness, $\mathregular{sone_{SHM}}$"))
            ax2.grid(alpha=0.075, linestyle='--')
            ax2.legend(bbox_to_anchor=(1, 0.85), title="Overall")

            # Filter signal to determine A-weighted time-averaged level
            if chan == 3:
                pA = A_weight_T(p_re, fs=sampleRate48k)
                LAeq2 = 20*np.log10(np.array([rms(pA[:, 0]),
                                              rms(pA[:, 1])])/2e-5)
                # take the higher channel level as representative (PD ISO/TS
                # 12913-3:2019 Annex D)
                LAeq = np.max(LAeq2)
                lr = np.argmax(LAeq2)
                # identify which channel is higher
                if lr == 0:
                    whichEar = " left ear"
                else:
                    whichEar = " right ear"
                # end of if branch to identify which channel is higher

                chan_lab = chan_lab + whichEar
            else:
                pA = A_weight_T(p_re[:, chan], fs=sampleRate48k)
                LAeq = 20*np.log10(rms(pA)/2e-5)

            fig.suptitle(t=(chan_lab + " signal sound pressure level = " +
                            str(roundTrad(LAeq, 1)) +
                            r"dB $\mathregular{\mathit{L}_{Aeq}}$"))
            fig.show()
        # end of for loop over channels
    # end of if branch for plotting

    # %% Output assignment

    # Assign outputs to structure
    if chansOut == 3:
        loudnessSHM = dict()
        loudnessSHM.update({'specLoudness': specLoudness[:, :, 0:2]})
        loudnessSHM.update({'specTonalLoudness': specTonalLoudness[:, :, 0:2]})
        loudnessSHM.update({'specNoiseLoudness': specNoiseLoudness[:, :, 0:2]})
        loudnessSHM.update({'specLoudnessPowAvg': specLoudnessPowAvg[:, 0:2]})
        loudnessSHM.update({'loudnessTDep': loudnessTDep[:, 0:2]})
        loudnessSHM.update({'loudnessPowAvg': loudnessPowAvg[0:2]})
        loudnessSHM.update({'specLoudnessBin': specLoudness[:, :, 3]})
        loudnessSHM.update({'specTonalLoudnessBin': specTonalLoudness[:, :, 3]})
        loudnessSHM.update({'specNoiseLoudnessBin': specNoiseLoudness[:, :, 3]})
        loudnessSHM.update({'specLoudnessPowAvgBin': specLoudnessPowAvg[:, 3]})
        loudnessSHM.update({'loudnessTDepBin': loudnessTDep[:, 3]})
        loudnessSHM.update({'loudnessPowAvgBin': loudnessPowAvg[3]})
        loudnessSHM.update({'bandCentreFreqs': bandCentreFreqs})
        loudnessSHM.update({'timeOut': timeOut})
        loudnessSHM.update({'soundField': soundField})
    else:
        loudnessSHM.update({'specLoudness': specLoudness})
        loudnessSHM.update({'specTonalLoudness': specTonalLoudness})
        loudnessSHM.update({'specNoiseLoudness': specNoiseLoudness})
        loudnessSHM.update({'specLoudnessPowAvg': specLoudnessPowAvg})
        loudnessSHM.update({'loudnessTDep': loudnessTDep})
        loudnessSHM.update({'loudnessPowAvg': loudnessPowAvg})
        loudnessSHM.update({'bandCentreFreqs': bandCentreFreqs})
        loudnessSHM.update({'timeOut': timeOut})
        loudnessSHM.update({'soundField': soundField})

# end of acousticSHLoudness function
