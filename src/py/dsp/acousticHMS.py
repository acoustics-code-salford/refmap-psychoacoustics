# -*- coding: utf-8 -*-
"""
acousticHMS.py
------------

Acoustic signal analysis routines implementing the Hearing Model of Sottek, as
defined in the ECMA-418-2 standard (currently 2022).

Requirements
------------
numpy (1.23.4)
scipy (1.9.3)
matplotlib (3.6.2)

Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford
 
Date created: 27/10/2023
Date last modified: 29/10/2023
Python version: 3.10.11

Copyright statement: This file and code is part of work undertaken within
the RefMap project (www.refmap.eu), and is subject to licence as detailed
in the code repository
(https://github.com/acoustics-code-salford/refmap-psychoacoustics)

Checked by:
Date last checked:

"""

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from scipy.signal import (bilinear, lfilter, lfilter_zi,
                          resample_poly, sosfilt, sosfreqz)
from math import gcd
from dsp.filterFuncs import A_weight_T

# set plot parameters
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams.update({'font.size': 14, 'mathtext.fontset': 'stix'})
mpl.rcParams['figure.autolayout'] = True


def acousticHMSPreProc(signal, blockSize, hopSize):
    """Function signalFadePad = acousticHMSPreProc(signal, blockSize, hopSize)

    Returns signal with fade-in and zero-padding pre-processing according to
    ECMA-418-2:2022 (the Hearing Model of Sottek) for an input signal.

    Inputs
    ------
    signal : 1D or 2D array
             the input signal/s

    blockSize : integer
                the maximum signal segmentation block size

    hopSize : integer
              the maximum signal segmentation hop size
              = (1 - overlap)*blockSize

    Returns
    -------
    signalFadePad : 1D or 2D array
                    the output faded, padded signal

    Assumptions
    -----------
    The input signal is oriented with time on axis 0 (and channel # on axis
    1), ie, the fade and padding operation is applied along axis 0.
    The input signal must be sampled at 48 kHz.

    Requirements
    ------------
    numpy

    Ownership and Quality Assurance
    -------------------------------
    Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
    Institution: University of Salford

    Date created: 26/09/2023
    Date last modified: 19/10/2023
    Python version: 3.10.11

    Copyright statement: This file and code is part of work undertaken within
    the RefMap project (www.refmap.eu), and is subject to licence as detailed
    in the code repository
    (https://github.com/acoustics-code-salford/refmap-psychoacoustics)

    Checked by:
    Date last checked:

    """
    # Signal processing

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

    # Fade in weighting function Section 5.1.2 ECMA-418-2:2022

    fadeWeight = 0.5 - 0.5*np.cos(np.pi*np.arange(0, 240)/240)[:, np.newaxis]
    # Apply fade in
    signalFade = np.concatenate((fadeWeight*signal[0:240, :],
                                 signal[240:, :]))

    # Zero-padding Section 5.1.2 ECMA-418-2:2022
    n_zeross = int(blockSize)  # start zero-padding
    n_samples = signal.shape[0]
    n_new = int(hopSize*(np.ceil((n_samples + hopSize + n_zeross)/hopSize) - 1))
    n_zerose = n_new - n_samples  # end zero-padding
    # Apply zero-padding
    signalFadePad = np.concatenate((np.zeros((n_zeross, numChans)), signalFade,
                                    np.zeros((n_zerose, numChans))))

    return signalFadePad  # end of acousticHMSPreProc function


def acousticHMSOutMidEarFilter(signal, outplot=False):
    """signalFiltered  = acousticHMSOutMidEarFilter_(signal, outplot)

    Returns signal filtered for outer and middle ear response according to
    ECMA-418-2:2022 (the Hearing Model of Sottek) for an input calibrated
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

    Requirements
    ------------
    numpy
    scipy

    Ownership and Quality Assurance
    -------------------------------
    Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
    Institution: University of Salford

    Date created: 29/10/2023
    Date last modified: 29/10/2023
    Python version: 3.10.11

    Copyright statement: This file and code is part of work undertaken within
    the RefMap project (www.refmap.eu), and is subject to licence as detailed
    in the code repository
    (https://github.com/acoustics-code-salford/refmap-psychoacoustics)

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

    return signalFiltered  # end of acousticHMSOutMidEarFilter function
