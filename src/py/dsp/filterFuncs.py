# -*- coding: utf-8 -*-
"""
filterFuncs.py
------------

Filter functions:

- Frequency weightings in the time domain
- Time weighting

Requirements
------------
numpy (1.26.3)
scipy (1.11.4)

Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford

Date created: 22/01/2024
Date last modified: 02/06/2025
Python version: 3.11.5

Copyright statements: This file and code is part of work undertaken within
the RefMap project (www.refmap.eu), and is subject to licence as detailed
in the code repository
(https://github.com/acoustics-code-salford/refmap-psychoacoustics)

Checked by:
Date last checked:

"""

import numpy as np
from scipy.signal import (bilinear, butter, freqz, lfilter, lfilter_zi,
                          resample_poly, sosfilt, sosfreqz)
from src.py.dsp.noct import noctf
from math import gcd


def A_weight_T(x, fs, axis=0, check=False):
    """
    Return time-domain-filtered signal according to standard sound frequency
    weighting 'A'.

    Implements IIR filter via bilinear transform. Includes pre-warping of
    analogue design frequencies to compensate for bilinear transform frequency
    distortion.

    Upsamples signals to 36kHz (if necessary)
    before processing, to ensure compliance with IEC 61672-1 class 1 acceptance
    limits.
    
    Resampling frequency and pre-warping defined according to [1].

    Parameters
    ----------
    x : 1D or 2D array
        contains the time signals to be weighted (filtered)
    fs : number
         the sampling frequency of the signals to be processed
    axis : integer
           the signal array axis along which to apply the filter
    check : boolean
            flag to check the filter against IEC acceptance limits

    Returns
    -------
    y : 1D or 2D array
        contains the weighted (filtered) time signals
    f : 1D array
        contains the frequencies (Hz) of the filter frequency response function
    H : 1D array
        contains the complex frequency response function values for each f

    Requirements
    ------------
    numpy
    scipy

    Assumptions
    -----------

    References
    ----------
    [1] Rimell, AN et al, 2015 - Design of digital filters for frequency
        weightings (A and C) required for risk assessments of workers exposed
        to noise. Industrial Health, 53, 21-27.

    """

    if fs < 36000:
        # upsampled sampling frequency
        fsu = 36000
        up = int(fsu/gcd(fsu, fs))
        down = int(fs/gcd(fsu, fs))
    else:
        fsu = fs
    dtu = 1/fsu

    G_Aw = 10**(2/20)
    
    f1 = np.sqrt((-1/(1 - np.sqrt(0.5))*(1e3**2 + 1e3*10**7.8/1e3**2
                                         - np.sqrt(0.5)*(1e3 + 10**7.8))
                  - np.sqrt((1/(1 - np.sqrt(0.5))*(1e3**2
                                                   + 1e3*10**7.8/1e3**2
                                                   - np.sqrt(0.5)*(1e3 + 10**7.8)))**2
                            - 4*1e3*10**7.8)) / 2)
    
    f4 = np.sqrt((-1/(1 - np.sqrt(0.5))*(1e3**2 + 1e3*10**7.8/1e3**2
                                         - np.sqrt(0.5)*(1e3 + 10**7.8))
                  + np.sqrt((1/(1 - np.sqrt(0.5))*(1e3**2
                                                   + 1e3*10**7.8/1e3**2
                                                   - np.sqrt(0.5)*(1e3 + 10**7.8)))**2
                            - 4*1e3*10**7.8)) / 2)
    
    f2 = 10**2.45*((3 - np.sqrt(5))/2)
    f3 = 10**2.45*((3 + np.sqrt(5))/2)

    w1 = 2*np.pi*f1
    w1w = 2/dtu*np.tan(w1*dtu/2)  # pre-warped frequency
    w4 = 2*np.pi*f4
    w4w = 2/dtu*np.tan(w4*dtu/2)  # pre-warped frequency
    w2 = 2*np.pi*f2
    w2w = 2/dtu*np.tan(w2*dtu/2)  # pre-warped frequency
    w3 = 2*np.pi*f3
    w3w = 2/dtu*np.tan(w3*dtu/2)  # pre-warped frequency

    B = np.array([G_Aw*w4w**2, 0, 0, 0, 0])
    A1 = [1.0, 2*w4w, (w4w)**2]
    A2 = [1.0, 2*w1w, (w1w)**2]
    A3 = [1.0, w3w]
    A4 = [1.0, w2w]
    A = np.convolve(np.convolve(np.convolve(A1, A2), A3), A4)

    b, a = bilinear(B, A, fsu)

    # determine filter initial conditions
    if len(x.shape) == 1 or axis == 1:
        zi = lfilter_zi(b, a)
    else:
        zi = lfilter_zi(b, a)[:, None]

    if check:
        # Check filter spec against acceptance limits
        f, H = freqz(b, a, worN=15000, fs=fsu)
        f = f[1:]  # discards zero-frequency
        H = H[1:]  # discards zero-frequency

        # IEC 61672-1:2013 acceptance limits (class 1)
        fm, _, _ = noctf(10, 20000, 3)
        Lm = np.array([-70.4, -63.4, -56.7, -50.5, -44.7, -39.4, -34.6, -30.2,
                      -26.2, -22.5, -19.1, -16.1, -13.4, -10.9, -8.6, -6.6, -4.8,
                      -3.2, -1.9, -0.8, 0.0, 0.6, 1.0, 1.2, 1.3, 1.2, 1.0, 0.5,
                      -0.1, -1.1, -2.5, -4.3, -6.6, -9.3])
        Ll = np.array([-9999.0, -9999.0, -4.0, -2.0, -1.5, -1.5, -1, -1, -1, -1,
                      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -0.7, -1, -1, -1,
                      -1, -1, -1, -1.5, -2, -2.5, -3, -5, -16, -9999]) + Lm
        Lu = np.array([3.0, 2.5, 2.0, 2.0, 2.0, 1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                      1, 1, 1, 1, 0.7, 1, 1, 1, 1, 1, 1, 1.5, 1.5, 1.5, 2, 2, 2.5,
                      3]) + Lm
        Lmf = np.interp(f, fm, Lm)  # Interpolate Lm onto f axis
        Luf = np.interp(f, fm, Lu)  # Interpolate Lu onto f axis
        Llf = np.interp(f, fm, Ll)  # Interpolate Lu onto f axis

    # Filter data on upsampled version

    if len(x.shape) == 1:
        # upsample signal (if necessary)
        if fsu > fs:
            x = resample_poly(x, up, down, padtype='line')
        # filter signal
        y, _ = lfilter(b, a, x, zi=zi*x[0])
        # if upsampled, downsample to original fs
        if fsu > fs:
            y = resample_poly(y, down, up, padtype='line')

    elif len(x.shape) == 2:
        # upsample signal (if necessary)
        if fsu > fs:
            x = resample_poly(x, up, down, axis=axis, padtype='line')
        # filter signal
        y, _ = lfilter(b, a, x, axis=axis, zi=zi*np.take(x, [0], axis=axis))
        # if upsampled, downsample to original fs
        if fsu > fs:
            y = resample_poly(y, down, up, axis, padtype='line')

    else:
        raise TypeError("\nInput must be 1d or 2d array")

    return y


def hrv_weight_T(x, fs, weight, axis=0, f_warp=True):
    """
    Return time-domain-filtered signal according to standard human response to
    vibration frequency weightings defined in ISO 8041:2005 (currently
    implemented: 'Wb', 'Wd').

    Implements IIR filter via bilinear transform. Includes pre-warping of
    analogue design frequencies to compensate for bilinear transform frequency
    distortion.

    Upsamples signals to minimum sample frequency (if necessary) according to
    [1] before processing, to ensure compliance with IEC 61672-1 class 1
    acceptance limits.

    Parameters
    ----------
    x : 1D or 2D array
        contains the time signals to be weighted (filtered)
    fs : number
        the sampling frequency of the signals to be processed
    weight : string
        indicates the weighting type as defined in BS EN ISO 8041:2005
    axis : integer
        the signal array axis along which to apply the filter
    f_warp : Boolean
        indicates whether frequency pre-warping should be applied to compensate
        filter coefficients for frequency distortion introduced by the bilinear
        transform

    Returns
    -------
    y : 1D or 2D array
        contains the weighted (filtered) time signals
    f : 1D array
        contains the frequencies (Hz) of the filter frequency response function
    H : 1D array
        contains the complex frequency response function values for each ff

    Requirements
    ------------
    numpy
    scipy

    Assumptions
    -----------

    References
    ----------
    [1] Rimell, AN et al, 2007 - Design of digital filters for frequency
        weightings required for risk assessments of workers exposed to
        vibration. Industrial Health, 45, 512-519.

    """
    weights = ['Wb', 'Wd']
    if weight.lower() not in [s.lower() for s in weights]:
        raise ValueError("\nWeighting must be selected from 'Wb' or 'Wd'")

    if weight.lower() == 'wb':
        # filter frequencies
        f1 = 0.4
        f2 = 100.0
        f3 = 16.0
        f4 = 16.0
        f5 = 2.5
        f6 = 4.0
        # filter resonant quality factors
        Q4 = 0.55
        Q5 = 0.9
        Q6 = 0.95
        # filter gain for zpk format (not currently used)
#        K = 1.024
        # check fs is above minimum required - if not, scale up fs for
        # resampling
        if fs < 12*f2:
            fsu = 12*f2
            up = int(fsu/gcd(fsu, fs))
            down = int(fs/gcd(fsu, fs))
        else:
            fsu = fs

    elif weight.lower() == 'wd':
        # filter frequencies
        f1 = 0.4
        f2 = 100.0
        f3 = 2.0
        f4 = 2.0
        f5 = np.inf
        f6 = np.inf
        # filter resonant quality factors
        Q4 = 0.63
        Q5 = 1.0
        Q6 = 1.0
        # filter gain for zpk format (not currently used)
#        K = 1.0
        # check fs is above minimum required - if not, scale up fs for
        # resampling
        if fs < 12*f2:
            fsu = 12*f2
            up = int(fsu/fs)
        else:
            fsu = fs

    # Nyquist frequency
#    fnyq = fsu/2
    # discrete time step from fs
    dtu = 1/fsu
    # if necessary, resample signals
    if fs < fsu:
        x = resample_poly(x, up, down, axis=axis, padtype='line')
    # convert filter frequencies to angular
    w3 = 2*np.pi*f3
    w4 = 2*np.pi*f4
    w5 = 2*np.pi*f5
    w6 = 2*np.pi*f6
    # apply frequency pre-warping to compensate for bilinear distortion
    if f_warp:
        if not np.isinf(w3):
            w3 = 2/dtu*np.tan(w3*dtu/2)
        if not np.isinf(w4):
            w4 = 2/dtu*np.tan(w4*dtu/2)
        if not np.isinf(w5):
            w5 = 2/dtu*np.tan(w5*dtu/2)
        if not np.isinf(w6):
            w6 = 2/dtu*np.tan(w6*dtu/2)

    # parameters for band-limiting filter
    sos1 = butter(2, f1, btype='highpass', output='sos', fs=fsu)
    sos2 = butter(2, f2, output='sos', fs=fsu)
    # band-limiting filter frequency response functions
    fN = 150000
    f, H1 = sosfreqz(sos1, worN=fN, fs=fsu)
    _, H2 = sosfreqz(sos2, worN=fN, fs=fsu)

    # parameters for a-v transition filter
    if not np.logical_and(f3 == np.inf, f4 == np.inf):
        B3 = np.array([1/w3, 1])
        A3 = np.array([1/w4**2, 1/Q4/w4, 1])
        b3, a3 = bilinear(B3, A3, fs=fsu)
        # a-v transition filter frequency response function
        _, H3 = freqz(b3, a3, worN=fN, fs=fsu)
    else:
        # ensure that infinite frequencies are handled
        H3 = 1

    # parameters for upward step filter
    if not np.logical_and(f5 == np.inf, f6 == np.inf):
        B4 = np.array([1/w5**2, 1/Q5/w5, 1])*w5**2/w6**2
        A4 = np.array([1/w6**2, 1/Q6/w6, 1])
        b4, a4 = bilinear(B4, A4, fs=fsu)
        # upward step filter frequency response function
        _, H4 = freqz(b4, a4, worN=fN, fs=fsu)
    else:
        # ensure that infinite frequencies are handled
        H4 = 1

    # total band limiting filter frequency response function
    Hbl = H1*H2
    # total frequency response function
    H = Hbl*H3*H4

    # apply filters
    y = sosfilt(sos2, x, axis=axis)  # lowpass
    y = sosfilt(sos1, y, axis=axis)  # highpass

    if not np.logical_and(f3 == np.inf, f4 == np.inf):
        y = lfilter(b3, a3, y, axis=axis)  # a-v transition

    if not np.logical_and(f5 == np.inf, f6 == np.inf):
        y = lfilter(b4, a4, y, axis=axis)  # upward step

    # remove last sample
    if axis == 0:
        y = y[:-1, :]
    elif axis == 1:
        y = y[:, :-1]

    # if necessary, resample signals back to original fs
    if fsu > fs:
        y = resample_poly(y, down, up, axis=axis, padtype='line')

    return y, f, H


def time_weight(x, fs, tau=0.125, axis=0):
    """
    Return signal filtered to give exponential time-weighted version of
    input signal x. The low-pass filter represents the root-mean-square
    integration function with exponential time-weighting. Common values for the
    input time constant tau are 0.125 ('Fast' weighting) and 1.0 ('Slow'
    weighting).

    Parameters
    ----------
    x : 1D or 2D array
        contains the time signals to be weighted (filtered)
    fs : number
         the sampling frequency of the signals to be processed
    tau : number, optional
          exponential time-weighting constant. The default is 0.125.
    axis : integer, optional
           the signal array axis along which to apply the filter. The default
           is 0.

    Returns
    -------
    y : 1D or 2D array
        contains the filtered output signal/s

    """

    wc = 1/tau
    dt = 1/fs
    k = np.exp(-wc*dt)
    # filter coefficients
    b = np.array([1 - k])
    a = np.array([1, -k])

    # generate filter with initial condition
    zi = lfilter_zi(b, a)
    
    # if necessary, expand zi dimension to avoid lfilter error
    if zi.ndim != x.ndim:
        zi = np.expand_dims(zi, axis)

    y, zf = lfilter(b, a, x**2, axis=axis, zi=zi*x[0]**2)  # apply filter
    y = np.sqrt(y)
    return y
