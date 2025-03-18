# -*- coding: utf-8 -*-
# %% Preamble
"""
acousticSHMSubs.py
--------------

Acoustic signal analysis subroutines for implementing the Sottek Hearing Model,
as defined in the ECMA-418-2 standard (currently 2024).

Requirements
------------
numpy
scipy
matplotlib

Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford
 
Date created: 27/10/2023
Date last modified: 04/02/2025
Python version: 3.11

Copyright statement: This file and code is part of work undertaken within
the RefMap project (www.refmap.eu), and is subject to licence as detailed
in the code repository
(https://github.com/acoustics-code-salford/refmap-psychoacoustics)

As per the licensing information, please be aware that this code is WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

Parts of this code were developed from an original MATLAB file
'SottekTonality.m' authored by Matt Torjussen (14/02/2022), based on
implementing ECMA-418-2:2020. The original code has been reused and translated
here with permission.

Checked by:
Date last checked:

"""

# %% Imports and parameter setup
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from scipy.signal import (bilinear, freqz, lfilter, lfilter_zi,
                          resample_poly, sosfilt, sosfreqz)
from scipy.special import (comb)
from math import gcd
from dsp.filterFuncs import A_weight_T

# set plot parameters
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams.update({'font.size': 14, 'mathtext.fontset': 'stix'})
mpl.rcParams['figure.autolayout'] = True


# %% shmAuditoryFiltBank
def shmAuditoryFiltBank(signal, outplot=False):
    """Function signalFiltered = shmAuditoryFiltBank(signal, outplot)

    Returns a set of signals, bandpass filtered for the inner ear response
    in each half-Bark critical band rate scale width, according to
    ECMA-418-2:2024 (the Sottek Hearing Model) for an input calibrated
    single (sound pressure) time-series signal.

    Inputs
    ------
    signal : 1D or 2D array
             the input signal as single audio (sound pressure) signal

    outplot : Boolean true/false (default: false)
              flag indicating whether to generate a figure a frequency and phase
              response figure for the filter bank

    Returns
    -------
    signalFiltered : 2D array
                     the filtered signals 

    Assumptions
    -----------
    The input signal is oriented with time on axis 1, ie, the filter
    operation is applied along axis 1.
    The input signal must be sampled at 48 kHz.

    Checked by:
    Date last checked:
    %
    # %Arguments validation
    #     arguments (Input)
    #         signal (:, 1) double {mustBeReal}
    #         outplot {mustBeNumericOrLogical} = false
    #     end
    """

    # %% Arguments validation

    # TODO - see above

    # %% Define constants

    sampleRate48k = 48e3  # Signal sample rate prescribed to be 48kHz (to be used for resampling), Section 5.1.1 ECMA-418-2:2024
    deltaFreq0 = 81.9289  # defined in Section 5.1.4.1 ECMA-418-2:2024
    c = 0.1618  # Half-Bark band centre-frequency demoninator constant defined in Section 5.1.4.1 ECMA-418-2:2024

    halfBark = np.arange(0.5, 27, 0.5)  # half-critical band rate scale
    bandCentreFreqs = (deltaFreq0/c)*np.sinh(c*halfBark)  # Section 5.1.4.1 Equation 9 ECMA-418-2:2024
    dfz = np.sqrt(deltaFreq0**2 + (c*bandCentreFreqs)**2)  # Section 5.1.4.1 Equation 10 ECMA-418-2:2024

    # %% Signal processing

    # Apply auditory filter bank
    # --------------------------
    # Filter equalised signal using 53 1/2Bark ERB filters according to
    # Section 5.1.4.2 ECMA-418-2:2024

    k = 5  # filter order = 5, footnote 5 ECMA-418-2:2024
    e_i = [0, 1, 11, 11, 1]  # filter coefficients for Section 5.1.4.2 Equation 15 ECMA-418-2:2024

    signalFiltered = np.zeros((len(signal), len(halfBark)))
    for zBand in range(53):
        # Section 5.1.4.1 Equation 8 ECMA-418-2:2024
        tau = (1/(2**(2*k - 1)))*comb(2*k - 2, k - 1)*(1/dfz[zBand])

        d = np.exp(-1./(sampleRate48k*tau))  # Section 5.1.4.1 ECMA-418-2:2024

        # Band-pass modifier Section 5.1.4.2 Equation 16/17 ECMA-418-2:2024
        bp = np.exp((1j*2*np.pi*bandCentreFreqs[zBand]*np.arange(0, k + 2))/sampleRate48k)

        # Feed-backward coefficients, Section 5.1.4.2 Equation 14 ECMA-418-2:2024
        m = range(1, k + 1)
        a_m = np.append(1, ((-d)**m)*comb(k, range(1, 6)))*bp[0:k + 1]

        # Feed-forward coefficients, Section 5.1.4.2 Equation 15 ECMA-418-2:2024
        m = range(0, k)
        i = range(1, k)
        b_m = ((((1 - d)**k)/np.sum(e_i[1:]*(d**i)))*(d**m)*e_i)*bp[0:k]

        # Recursive filter Section 5.1.4.2 Equation 13 ECMA-418-2:2024
        # Note, the results are complex so 2x the real-valued band-pass signal
        # is required.
        signalFiltered[:, zBand] = 2*np.real(lfilter(b_m, a_m, signal))

        # Plot figures

        if outplot:
            f, H = freqz(b_m, a=a_m, worN=int(10e3), whole=True, fs=48e3)
            phir = np.angle(H)
            phirUnwrap = np.unwrap(p=phir, discont=np.pi, axis=0)
            phiUnwrap = phirUnwrap/np.pi*180
            # Plot frequency and phase response for filter
            if zBand == 0:
                fig, axs = plt.subplots(nrows=2, ncols=1)
                ax1 = axs[0]
                ax2 = axs[1]
                ax1.semilogx(f, 20*np.log10(abs(H)))
                ax2.semilogx(f, phiUnwrap)

            if zBand == 52:
                ax1.set(xlim=[20, 20e3],
                        xticks=[31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],
                        xticklabels=["31.5", "63", "125", "250", "500", "1k",
                                     "2k", "4k", "8k", "16k"],
                        xlabel="Frequency, Hz", ylabel="$H$, dB")
                ax1.minorticks_off()
                ax1.grid(alpha=0.15, linestyle='--')

                ax2.set(xlim=[20, 20e3],
                        xticks=[31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],
                        xticklabels=["31.5", "63", "125", "250", "500", "1k",
                                     "2k", "4k", "8k", "16k"],
                        xlabel="Frequency, Hz", ylabel=r"Phase angle$^\degree$")
                ax2.minorticks_off()
                ax2.grid(alpha=0.15, linestyle='--')

    return signalFiltered  # end of shmAuditoryFiltBank function


# %% ShmBasisLoudness
def ShmBasisLoudness(signalSegmented, bandCentreFreq):
    """Function ShmBandBasisLoudness(signalSegmented, bandCentreFreq)
    
    Returns rectified input and basis loudness in specified half-Bark
    critical band according to ECMA-418-2:2024 (the Sottek Hearing Model)
    for an input band-limited signal, segmented into processing blocks
    
    Inputs
    ------
    signalSegmented : 2D or 3D matrix
              input band-limited segmented signal(s)
    
    bandCentreFreq : double (optional, default = [])
             half-Bark critical band centre frequency - if empty, all
             bands are assumed to be present in the input segmented
             signal matrix
    
    Returns
    -------
    signalRectSeg : 2D or 3D matrix
            rectified band-limited segmented signal, orientated as
            per the input
    
    basisLoudness : 2D or 3D matrix
            basis loudness in each block, orientated as per the input
    
    blockRMS : column vector or 2D matrix
       RMS for each block, orientated as per the input with singleton
       dimension removed
    
    Assumptions
    -----------
    The input signal is a segmented signal (either band-limited, or arranged
    with half-Bark critical bands over the third dimension) obtained using
    acousticSHMAuditoryFiltBank.m and ShmSignalSegment.m
    
    Requirements
    ------------
    None
    
    Ownership and Quality Assurance
    -------------------------------
    Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk) &
     Matt Torjussen (matt@anv.co.uk)
    Institution: University of Salford / ANV Measurement Systems
    
    Date created: 27/09/2023
    Date last modified: 21/10/2024
    MATLAB version: 2023b
    
    Copyright statement: This file and code is part of work undertaken within
    the RefMap project (www.refmap.eu), and is subject to licence as detailed
    in the code repository
    (https://github.com/acoustics-code-salford/refmap-psychoacoustics)
    
    As per the licensing information, please be aware that this code is
    WITHOUT ANY WARRANTY without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    
    This code was developed from an original file 'SottekTonality.m' authored
    by Matt Torjussen (14/02/2022), based on implementing ECMA-418-2:2020.
    The original code has been reused and updated here with permission.
    
    Checked by:
    Date last checked:
    """

    # %% Arguments validation
        # arguments (Input)
        #     signalSegmented double {mustBeReal}
        #     bandCentreFreq double {mustBePositive} = []
        # end
    # TODO
    
    # check if input is 2D and includes band centre frequency - otherwise raise
    # error
    if isempty(bandCentreFreq) && length(size(signalSegmented)) == 2
        raise ValueError("Band centre frequency must be specified for single band-limited input signal")
    end
    
    # check if input band centre frequency is not a vector (arguments
    # validation does not allow empty default with specified size)
    if ~isempty(bandCentreFreq) && max(size(bandCentreFreq)) ~= 1
        error("Band centre frequency input must be a single value")
    end
    
    # %% Define constants
    
    deltaFreq0 = 81.9289  # defined in Section 5.1.4.1 ECMA-418-2:2024
    c = 0.1618  # Half-Bark band centre-frequency denominator constant defined in Section 5.1.4.1 ECMA-418-2:2024
    
    halfBark = 0.5:0.5:26.5  # half-critical band rate scale
    bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark)  # Section 5.1.4.1 Equation 9 ECMA-418-2:2024
    
    cal_N = 0.0211668  # Calibration factor from Section 5.1.8 Equation 23 ECMA-418-2:2024
    cal_Nx = 1.00132  # Calibration multiplier (Footnote 8 ECMA-418-2:2024)
    #cal_N*cal_Nx = 0.021194740176  # Adjusted calibration factor
    #cal_Nx = 1.001398416387928  # Calibration multiplier (Footnote 8 ECMA-418-2:2024)
    #cal_N*cal_Nx = 0.0211964  # Adjusted calibration factor
    
    a = 1.5  # Constant (alpha) from Section 5.1.8 Equation 23 ECMA-418-2:2024
    
    # Values from Section 5.1.8 Table 2 ECMA-418-2:2024
    p_threshold = 2e-5*10.^((15:10:85)/20).T
    v = [1, 0.6602, 0.0864, 0.6384, 0.0328, 0.4068, 0.2082, 0.3994, 0.6434]
    
    # Loudness threshold in quiet Section 5.1.9 Table 3 ECMA-418-2:2024
    LTQz = [0.3310, 0.1625, 0.1051, 0.0757, 0.0576, 0.0453, 0.0365, 0.0298,...
            0.0247, 0.0207, 0.0176, 0.0151, 0.0131, 0.0115, 0.0103, 0.0093,...
            0.0086, 0.0081, 0.0077, 0.0074, 0.0073, 0.0072, 0.0071, 0.0072,...
            0.0073, 0.0074, 0.0076, 0.0079, 0.0082, 0.0086, 0.0092, 0.0100,...
            0.0109, 0.0122, 0.0138, 0.0157, 0.0172, 0.0180, 0.0180, 0.0177,...
            0.0176, 0.0177, 0.0182, 0.0190, 0.0202, 0.0217, 0.0237, 0.0263,...
            0.0296, 0.0339, 0.0398, 0.0485, 0.0622]
    
    # %% Input check
    
    if ~isempty(bandCentreFreq) && ~ismember(bandCentreFreq, bandCentreFreqs)
        error("Input half-Bark critical rate scale band centre frequency does not match ECMA-418-2:2024 values")
    end
    
    # %% Signal processing
    
    # Half Wave Rectification
    # -----------------------
    # Section 5.1.6 Equation 21 ECMA-418-2:2020
    signalRectSeg = signalSegmented
    signalRectSeg(signalSegmented <= 0) = 0
    
    # Calculation of RMS
    # ------------------
    # Section 5.1.7 Equation 22 ECMA-418-2:2024
    blockRMS = sqrt((2/size(signalRectSeg, 1))*sum(signalRectSeg.^2, 1))
    
    # Transformation into Loudness
    # ----------------------------
    # Section 5.1.8 Equations 23 & 24 ECMA-418-2:2024
    bandLoudness = cal_N*cal_Nx*(blockRMS/20e-6).*prod((1 + (blockRMS./p_threshold).^a).^((diff(v)/a)'))
    
    # remove singleton dimension from block RMS output
    blockRMS = squeeze(blockRMS)
    
    # Section 5.1.9 Equation 25 ECMA-418-2:2024
    if ~isempty(bandCentreFreq) && length(size(signalSegmented)) == 2
        # half-Bark critical band basis loudness
        basisLoudness = bandLoudness - LTQz(bandCentreFreq == bandCentreFreqs)
        basisLoudness(basisLoudness < 0) = 0
    else
        # basis loudness for all bands
        basisLoudness = bandLoudness - repmat(reshape(LTQz, [1, 1, 53]),...
                                                  1, size(bandLoudness, 2), 1)
        basisLoudness(basisLoudness < 0) = 0
    end
    
    # end of function


# %% shmPreProc
def shmPreProc(signal, blockSize, hopSize, padStart, padEnd):
    """Function signalFadePad = shmPreProc(signal, blockSize, hopSize,
                                           padStart, padEnd)

    Returns signal with fade-in and zero-padding pre-processing according to
    ECMA-418-2:2024 (the Sottek Hearing Model) for an input signal.

    Inputs
    ------
    signal : 1D or 2D array
             the input signal/s

    blockSize : integer
                the maximum signal segmentation block size

    hopSize : integer
              the maximum signal segmentation hop size
              = (1 - overlap)*blockSize

    padStart : Boolean
               flag to indicate whether to pad the start of the signal

    padEnd : Boolean
             flag to indicate whether to pad the end of the signal
 

    Returns
    -------
    signalFadePad : 1D or 2D array
                    the output faded, padded signal

    Assumptions
    -----------
    The input signal is oriented with time on axis 0 (and channel # on axis
    1), ie, the fade and padding operation is applied along axis 0.
    The input signal must be sampled at 48 kHz.

    Copyright statement: This file and code is part of work undertaken within
    the RefMap project (www.refmap.eu), and is subject to licence as detailed
    in the code repository
    (https://github.com/acoustics-code-salford/refmap-psychoacoustics)

    Checked by:
    Date last checked:

    """
    # %% Arguments validation

    # arguments (Input)
    #     signal (:, :) double {mustBeReal}
    #     blockSize (1, 1) {mustBeInteger}
    #     hopSize (1, 1) {mustBeInteger}
    #     padStart {mustBeNumericOrLogical} = true
    #     padEnd {mustBeNumericOrLogical} = true
    # end

    # TODO

    # %% Signal processing

    # Input pre-processing
    # --------------------
    #
    # Check signal dimensions and add axis if 1D input
    #
    if np.size(signal.shape) == 1:
        numChans = 1
        signal = signal[:, np.newaxis]
    else:
        numChans = signal.shape(1)

    if numChans > 2:
        raise ValueError("Input signal must be 1- or 2-channel")

    # Fade in weighting function Section 5.1.2 ECMA-418-2:2024

    fadeWeight = 0.5 - 0.5*np.cos(np.pi*np.arange(0, 240)/240)[:, np.newaxis]
    # Apply fade in
    signalFade = np.concatenate((fadeWeight*signal[0:240, :],
                                 signal[240:, :]))

    # Zero-padding Section 5.1.2 ECMA-418-2:2024
    if padStart:
        n_zeross = int(blockSize)  # start zero-padding
    else:
        n_zeross = 0

    if padEnd:
        n_samples = signal.shape[0]
        n_new = int(hopSize*(np.ceil((n_samples
                                      + hopSize
                                      + n_zeross)/hopSize) - 1))
        n_zerose = n_new - n_samples  # end zero-padding

    # Apply zero-padding
    signalFadePad = np.concatenate([np.zeros((n_zeross, numChans)),
                                    signalFade,
                                    np.zeros((n_zerose, numChans))])

    return signalFadePad  # end of shmPreProc function


def shmOutMidEarFilter(signal, outplot=False):
    """signalFiltered  = shmOutMidEarFilter_(signal, outplot)

    Returns signal filtered for outer and middle ear response according to
    ECMA-418-2:2024 (the Sottek Hearing Model) for an input calibrated
    audio (sound pressure) time-series signal. Optional plot output shows
    frequency and phase response of the filter.

    Inputs
    ------
    signal : 1D or 2D array
             the input signal (sound pressure)

    outplot : Boolean true/false (default: false)
              flag indicating whether to generate a frequency and phase
              response figure for the filter

    Returns
    -------
    signalFiltered : 1D or 2D array
                     the output filtered signal

    Assumptions
    -----------
    The input signal is oriented with time on axis 0 (and channel # on axis
    1), ie, the filtering operation is applied along axis 0.
    The input signal must be sampled at 48 kHz.

    Checked by:
    Date last checked:

    %% Arguments validation
        arguments (Input)
            signal (:, :) double {mustBeReal}
            outplot {mustBeNumericOrLogical} = false
        end
    """
    # Arguments validation
    if not type(signal) is np.ndarray:
        try:
            np.array(signal)
        except Exception:  # TODO specific exception handling
            raise TypeError("\nInput signal does not appear to be an array")

    if not isinstance(outplot, bool):
        raise ValueError("\nInput argument 'outplot' must be logical True/False")

    # Signal processing

    # Apply outer & middle ear filter bank
    # ------------------------------------
    #
    # Filter coefficients from Section 5.1.3.2 Table 1 ECMA-418-2:2022
    # b_0k = [1.015896, 0.958943, 0.961372, 2.225804, 0.471735, 0.115267,
    #         0.988029, 1.952238]
    # b_1k = [-1.925299, -1.806088, -1.763632, -1.43465, -0.366092, 0.0,
    #         -1.912434, 0.16232]
    # b_2k = [0.922118, 0.876439, 0.821788, -0.498204, 0.244145, -0.115267,
    #         0.926132, -0.667994]
    # a_0k = ones(size(b_0k))
    # a_1k = [-1.925299, -1.806088, -1.763632, -1.43465, -0.366092, -1.796003,
    #         -1.912434, 0.16232]
    # a_2k = [0.938014, 0.835382, 0.78316, 0.727599, -0.28412, 0.805838,
    #         0.914161, 0.284244]
    #
    # Accurate coefficient values
    b_0k = [1.01589602025559, 0.958943219304445, 0.961371976333197,
            2.22580350360974, 0.471735128494163, 0.115267139824401,
            0.988029297230954, 1.95223768730136]
    b_1k = [-1.92529887777608, -1.80608801184949, -1.76363215433825,
            -1.43465048479216, -0.366091796830044, 0.0, -1.91243380293387,
            0.162319983017519]
    b_2k = [0.922118060364679, 0.876438777856084, 0.821787991845146,
            -0.498204282194628, 0.244144703885020, -0.115267139824401,
            0.926131550180785, -0.667994113035186]
    a_0k = np.ones(len(b_0k))
    a_1k = [-1.92529887777608, -1.80608801184949, -1.76363215433825,
            -1.43465048479216, -0.366091796830044, -1.79600256669201,
            -1.91243380293387, 0.162319983017519]
    a_2k = [0.938014080620272, 0.835381997160530, 0.783159968178343,
            0.727599221415107, -0.284120167620817, 0.805837815618546,
            0.914160847411739, 0.284243574266175]

    # form 2nd-order section for transfer function (copy is needed to enforce
    # C-contiguity in the tranposed array)
    sos = np.array([b_0k, b_1k, b_2k, a_0k, a_1k, a_2k]).T.copy(order='C')

    # Section 5.1.3.2 ECMA-418-2:2022 Outer and middle/inner ear signal filtering
    signalFiltered = sosfilt(sos, signal, axis=0)

    # Plot figures

    if outplot:
        f, H = sosfreqz(sos, worN=10000, whole=True, fs=48e3)
        phir = np.angle(H, deg=False)
        phirUnwrap = np.unwrap(phir, discont=np.pi)
        phiUnwrap = phirUnwrap/np.pi*180
        # Plot frequency and phase response for filter
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 7))
        ax1.semilogx(f[1:], 20*np.log10(np.abs(H[1:])), color=[0.0, 0.2, 0.8])
        ax1.set_xticks(ticks=[31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],
                       labels=["31.5", "63", "125", "250", "500", "1k", "2k",
                               "4k", "8k", "16k"])
        ax1.axis(xmin=20, xmax=20e3, ymin=-30, ymax=10)
        ax1.tick_params(axis='x',          # changes apply to the x-axis
                        which='minor',      # minor ticks are affected
                        bottom=False,     # ticks along the bottom edge are off
                        top=False,         # ticks along the top edge are off
                        labelbottom=False)  # labels along the bottom edge are off

        ax1.set_xlabel("Frequency, Hz", fontname='Arial', fontsize=12)
        ax1.set_ylabel("$H$, dB", fontname='Arial', fontsize=12)
        ax1.grid(True, which='major', linestyle='--', alpha=0.5)

        ax2.semilogx(f[1:], phiUnwrap[1:], color=[0.8, 0.1, 0.8])
        ax2.set_xticks(ticks=[31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],
                       labels=["31.5", "63", "125", "250", "500", "1k", "2k",
                               "4k", "8k", "16k"])
        ax2.set_yticks(ticks=np.arange(-180, 210, 30))
        ax2.axis(xmin=20, xmax=20e3, ymin=-180, ymax=90)
        ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
        ax2.tick_params(axis='x',          # changes apply to the x-axis
                        which='minor',      # minor ticks are affected
                        bottom=False,     # ticks along the bottom edge are off
                        top=False,         # ticks along the top edge are off
                        labelbottom=False)  # labels along the bottom edge are off
        ax2.set_xlabel("Frequency, Hz", fontname='Arial', fontsize=12)
        ax2.set_ylabel("Phase angle, degrees", fontname='Arial', fontsize=12)
        ax2.grid(True, which='major', linestyle='--', alpha=0.5)

    return signalFiltered  # end of shmOutMidEarFilter function
