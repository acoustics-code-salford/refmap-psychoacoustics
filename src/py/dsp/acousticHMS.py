# -*- coding: utf-8 -*-
"""
acousticHMS.py
------------

Module provides acoustic signal analysis routines
implementing the Hearing Model of Sottek, as defined
in the ECMA-418-2 standard (currently 2022).

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
Date last modified: 27/10/2023
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
    Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
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
