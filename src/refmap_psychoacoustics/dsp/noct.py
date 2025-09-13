# -*- coding: utf-8 -*-
"""
noct.py
------------

Fractional octave band filtering according to the class 1 specification of
IEC 61260:2014

Requirements
------------
numpy (1.23.4)
scipy (1.9.3)

Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford
 
Date created: 29/10/2023
Date last modified: 05/11/2024
Python version: 3.11

Copyright statement: This file and code is part of work undertaken within
the RefMap project (www.refmap.eu), and is subject to licence as detailed
in the code repository
(https://github.com/acoustics-code-salford/refmap-psychoacoustics)

Checked by:
Date last checked:

"""

from scipy.signal import (butter, resample_poly, sosfilt, sosfiltfilt, sosfilt_zi,
                          sosfreqz)
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Arial'


def noctf(fl, fh, n):
    """
    Return consecutive range of (exact) mid-frequencies, and lower and upper
    band-edge frequencies for 1/n fractional-octave-band filters according to
    BS EN 61260:2014

    Parameters
    ----------
    fl : number
      The lowest frequency band of interest
    fh : number
      The highest frequency band of interest
    n : integer
      A number defining the octave fraction 1/n

    Returns
    -------
    fm : numpy array
      A 1D array of the mid frequencies
    f1 : numpy array
      A 1D array of the lower band-edge frequencies
    f2 : numpy array
      A 1D array of the upper band-edge frequencies

    Assumptions
    -----------

    """
    if not isinstance(n, int):
        raise TypeError("\nOctave fraction denominator must be an integer")

    ind = np.arange(-15*n, 15*n, 1)  # range of frequency indices
    G10 = 10**(3/10)  # octave ratio coefficient (base-ten)
    OctRatio = G10**(0.5/n)  # octave ratio

    if np.mod(n, 1) == 0:
        f = G10**(ind/n)*1000
    else:
        f = G10**((2*ind + 1)/(2*n))*1000

    il = np.abs(f - fl).argmin()  # find nearest lower exact frequency
    ih = np.abs(f - fh).argmin()  # find nearest higher exact frequency

    fm = f[il:ih+1]  # output range of exact fractional frequencies
    f1 = fm/OctRatio  # output range of exact lower band-edge frequencies
    f2 = fm*OctRatio  # output range of exact upper band-edge frequencies

    return fm, f1, f2


def noctfnoct(fl, fh, n):
    """
    Return consecutive range of nominal mid-frequencies for 1/1 or 1/3
    fractional-octave-band filters according to BS EN 61260:2014

    Parameters
    ----------
    fl : number
      The lowest frequency band of interest
    fh : number
      The highest frequency band of interest
    n : integer
      A number defining the octave fraction 1/n (n = 1 or n = 3)

    Returns
    -------
    fnm : numpy array
      A 1D array of the nominal mid frequencies

    Assumptions
    -----------

    """
    if not isinstance(n, int):
        raise TypeError("\nOctave fraction denominator must be an integer")

    if n not in [1, 3]:
        raise ValueError("\nOctave fraction denominator must take a value of "
                         "1 or 3")

    ind = np.arange(-15*n, 15*n, 1)  # range of frequency indices
    G10 = 10**(3/10)  # octave ratio coefficient (base-ten)

    if np.mod(n, 1) == 0:
        f = G10**(ind/n)*1000
    else:
        f = G10**((2*ind + 1)/(2*n))*1000

    il = np.abs(f - fl).argmin()  # find nearest lower exact frequency
    ih = np.abs(f - fh).argmin()  # find nearest higher exact frequency

    fm = f[il:ih+1]  # output range of nominal fractional frequencies

    fn = np.array([0.0315, 0.04, 0.05, 0.063, 0.08, 0.1, 0.125, 0.16, 0.2,
                   0.25])
    
    fn = (np.append(np.append(np.append(np.append(np.append(
          np.append(np.append(np.append(fn, fn*10), fn*100), fn*1e3), fn*1e4),
          fn*1e5), fn*1e6), fn*1e7), fn*1e8))

    im = np.zeros(len(fm), dtype=int)
    for ii, ff in enumerate(fm):
        im[ii] = np.abs(ff - fn).argmin()

    fnm = fn[im]

    return fnm


def noctf_fbw(f, fbw_order, fs, btype):
    """
    Return filter cutoff frequencies adjusted (pre-warped) to compensate for
    zero-phase forwards-backwards Butterworth filter response
    
    Developed from filter response functions for (one-way) orders 1-64
    (ie two-way orders 2-128) using linear regression over values adjusted to
    give desired output

    Parameters
    ----------
    
    f : float or array of floats
        the input frequency which is to be adjusted/warped. If btype='bandpass'
        f must be list of two sets of frequencies, the lower band cutoff and
        upper band cutoff frequencies (can be a list of arrays)
    fbw_order : integer
        the design target order of the fowards-backwards filter, ie 2x the
        value of the order input to the Butterworth filter function
    fs : integer
         the sampling frequency
    btype : string
        the type of filter; choose from 'lowpass', 'highpass' or 'bandpass'
        
    Returns
    -------
    f : float or array of floats
        Adjusted frequencies. If 'btype' is 'highpass' or 'lowpass', f will be
        the value or array of values for the cutoff frequencies. If 'btype' is
        'bandpass', two sets of f (f1, f2) are returned corresponding with the
        upper and lower band-edge frequencies

    """
    if not isinstance(fbw_order, int):
        raise TypeError("\nForwards-backwards filter order parameter "
              "'fbw_order' must be an integer")

    if fbw_order < 2 or fbw_order > 128:
        raise ValueError("\nForwards-backwards filter order parameter "
              "'fbw_order' must take a value in the range 2-128")

    K = np.log(fbw_order)
    if btype == 'highpass':
        c = np.array([-1.95309206895014e-03, 3.3409655894888e-02,
                      -2.27454951482515e-01, 7.88665396386179e-01,
                      -1.46009481744288, 2.24142857142864])
        f_fbw = f/(np.sum(c*K**np.arange(5, -1, -1)))
        return f_fbw
    elif btype == 'lowpass':
        c = np.array([-1.22393769654216e-03, 2.17097125485945e-02,
                      -1.55752187929785e-01, 5.80441774578037e-01,
                      -1.17262034333309, 2.081])
        f_fbw = f*(np.sum(c*K**np.arange(5, -1, -1)))
        return f_fbw
    elif btype == 'bandpass':
        c = np.array([-3.5676481792823e-05, 1.02690386340958e-03,
                      -1.05748495945522E-02, 5.20134108353335E-02,
                      -1.27978679123245E-01, 1.13201])
        f1_fbw = f[0]/(np.sum(c*K**np.arange(5, -1, -1)))
        f2_fbw = f[1]*(np.sum(c*K**np.arange(5, -1, -1)))
        return f1_fbw, f2_fbw
    else:
        raise ValueError("\nArgument 'btype' must be 'highpass', 'lowpass' or 'bandpass'")


def noctlimits(n, fm):
    """
    Return limits for class 1 fractional octave band filter

    Calculates tolerance limits for fractional octave band filter to class 1
    as specified in BS EN IEC 61260:2014.
    
    Parameters
    ----------
    n : integer
        A number defining the octave fraction 1/n
    fm : sequence of floats
         Mid-frequencies of 1/n octave filters. Float, or list or numpy array
         (1D or singleton-2D) of floats

    Returns
    -------
    Lu : 1D array
        Upper limit level (in dB)
    Ll : 1D array
        Lower limit level (in dB)
    f : 1D or 2D array
        Frequency axis for limit levels as numpy array. If fm contains 
        several frequencies, each row in f (2D array) corresponds to a
        consecutive value in fm (1D array).


    """
    # Argument validation
    if n < 0 or not isinstance(n, int):
        raise ValueError("\nOctave fraction denominator must be a positive integer")

    if isinstance(fm, list):
        fm = np.array(fm)

    G = 10**(3/10)  # Octave frequency ratio (base-10)

    # Calculate normalised frequency break points
    fnorm = np.array([G**(1/8), G**(1/4), G**(3/8), G**(1/2)*(1 - 1e-16),
                          G**(1/2), G, G**2, G**3, G**4])
    fnorm = 1 + ((G**(1/2/n) - 1)/(np.sqrt(G) - 1))*(fnorm - 1)

    if np.size(fm) == 1:
        f = fm*np.append(np.append(np.flip(1/fnorm), 1), fnorm)
    else:
        f = np.zeros((19, len(fm)))
        for ii in range(0, len(fm)):
            f[ii, :] = fm[ii]*np.append(np.append(np.flip(1/fnorm), 1),
                                        fnorm)

    # Create limits
    Lu = np.array([-70, -60, -40.5, -16.6, -1.2, .4, .4, .4, .4, .4, .4,
                   .4, .4, .4, -1.2, -16.6, -40.5, -60, -70])
    Ll = np.array([-9999.0, -9999.0, -9999.0, -9999.0, -9999.0, -5.3, -1.4,
                   -.7, -.5, -.4, -.5, -.7, -1.4, -5.3, -9999.0, -9999.0,
                   -9999.0, -9999.0, -9999.0])

    return Lu, Ll, f


def noctfiltc(n, fm, fs, order=4, fwd_bwd=False, down_fact=1, plot=False):
    """
    Return second-order-sections (sos) for a 1/n fractional-octave-band filter
    according to BS EN IEC 61260:2014
        
    Parameters
    ----------
    n : integer
        A number defining the octave fraction 1/n
    fm : sequence of floats
         Mid-frequencies of 1/n octave filters. Float, or list or numpy array
         (1D or singleton-2D) of floats
    fs : integer
         the sampling frequency
    order : integer
            the filter order
    fwd_bwd : Boolean (default: False)
              determines whether the filter is designed as one-way or two-way
              (forwards-backwards)
    down_fact : integer
                downsampling factor to apply to filter design
                (TODO: currently unused)
    plot : Boolean
           determines whether a plot of the filter frequency response is
           displayed

    Returns
    -------
    sos : array
          the filter coefficients in second-order section form
    zi : float
         the initial condition for the filter

    """
    N = [1, 2, 3, 6, 12, 24, 48]
    if n not in N:
        raise ValueError(("\nOctave fraction denominator must take a value from: "
                          + str(N)))

    # Check and force order to be integer with warning
    if not isinstance(order, int):
        ordernonint = order
        order = int(order)
        print(("\nCaution: Butterworth filter must take integer order. Input order "
               + str(ordernonint) + " changed to integer order " + str(order)))

    if fwd_bwd:
        if np.remainder(order, 2) > 0:
            inorder = order
            order = order + 1
            print(("\nCaution: for forward-backward processing, even input " \
                   "(target) order is required. Odd input order "
                   + str(inorder) + " changed to even order " + str(order)))

    G10 = 10**(3/10)
    octRatio = G10**(0.5/n)

    f1 = fm/octRatio
    f2 = fm*octRatio

    # TODO: implement downsampling section that works with normal and fbw filtering
    # (will also need to be passed through to noct_filter)
    # # downsampling based on Nyquist and filter upper edge frequency
    # down_fact = np.maximum(np.minimum(np.floor(fs/(2 + 0.1))/f2, 50), 1).astype(int)

    # if down_fact > 1:
    #     fs_filt = fs/down_fact
    # else:
    #     fs_filt = fs

    if fwd_bwd is True:
        # Suggested pre-warp for 2nd order Butterworth filtfilt cutoff
        # from Robertson et al 2014 Research methods in Biomechanics (2e)
        # Seems to give erroneous results, even for 2nd order filter
        # passes = 2
        # C = ((2**(1/passes))-1)**(1/4)
        # f1 = f1*C
        # f2 = f2/C
        fbw_order = order
        order = int(order/2)

        # Pre-warping of cutoff frequencies to compensate for forward-backward
        f1, f2 = noctf_fbw([f1, f2], fbw_order, fs, 'bandpass')  # fs_filt, 'bandpass')

    sos = butter(order, [f1, f2], btype='bandpass', output='sos', fs=fs)  # fs=fs_filt)
    # initial condition for steady-state of step response to pass to filter function
    zi = sosfilt_zi(sos)

    if plot:
        ff, Hf = sosfreqz(sos, worN=150000, fs=fs)  # fs=fs_filt)  # Filter frequency response, 150000 values
        if fwd_bwd:
            Hf = Hf*np.conj(Hf)
        plt.plot(ff[1:], 20*np.log10(abs(Hf[1:])), label=str(np.rint(fm)))
        plt.xlim((fm/octRatio/3, fm*octRatio*3))
        plt.ylim((-100, 10))
        plt.grid(b=True, which='both')
        plt.xlabel("Frequency, Hz")
        plt.ylabel("Filter response, dB")
        plt.legend(loc=4)

    return sos, zi #, fs_filt  # TODO


def noct_filter(x, n, fm, fs, order=4, axis=0, fwd_bwd=False, check=True):
    """
    Return 1/n fractional-octave-band-filtered signals for a single passband
    using second-order-sections.

    Parameters
    ----------
    x : numpy array
        1D or 2D array of the input time signals to be filtered
    n : integer
        Number defining the octave fraction 1/n
    fm : number
        Mid-frequency for the fractional-octave-band filter
    fs : number
        Sampling frqeuency of the input signal
    order : integer
        Intended order for Butterworth IIR filter design
    axis : integer
        Axis of x along which to apply filter
    fwd_bwd : Boolean
        Indicates forward-backward (zero-phase) filtering to be used.
        In this case, 'order' is interpreted as the overall desired target
        filter order, and a caution is raised if 'order' is odd (as the target
        order must be even if forward-backward filtering is used). The filter
        is padded by a constant length taken from the smaller of: i) the length
        of x along 'axis' minus 1, or ii) a value calculated according to [2].
    check : Boolean
        Flag to check for filter compliance with BS EN IEC 61260:2014 [1].

    Returns
    -------
    y : numpy array
      An array of the filtered output signals with shape: x.shape

    Assumptions
    -----------
    x is a 2D array (default behaviour assumes axis for filter is 0)

    References
    ----------
    1. BSI, 2014. BS EN IEC 61260:2014
    2. Boore, DM, 2005. On pads and filters: processing strong-motion data. In:
        Bulletin of the Seismological Society of America, 95(2), 745-750.
    """
    N = [1, 2, 3, 6, 12, 24, 48]
    if n not in N:
        raise ValueError(("\nOctave fraction denominator must take a value from: "
                          + str(N)))

    fm, f1, f2 = noctf(fm, fm, n)  # convert to exact filter mid-frequency

    # derive filter coefficients
    sos, _ = noctfiltc(n, fm, fs, order=order, fwd_bwd=fwd_bwd)

    if check is True:
        # Check that filter corresponds to standard limits, or issue warning
        # The filter shape is checked in freq. range fc/5 < f < 5*fc
        Lu, Ll, f = noctlimits(n, fm)  # Calculate class 1 limits
        # Filter frequency response, 150000 values
        ff, Hf = sosfreqz(sos, worN=150000, fs=fs)

        if fwd_bwd:
            Hf = Hf*np.conj(Hf)

        Hf = Hf[1:]  # Remove zero frequency
        ff = ff[1:]

        Luff = np.interp(ff, f, Lu)  # Interpolate Lu onto ff axis
        # Interpolate Ll onto ff axis
        Llff = np.interp(ff, f, Ll, left=-999, right=-999)
        Llff[np.isnan(Llff)] = -999
        idx = np.nonzero(np.logical_and.reduce((ff > fm/5, ff < 5*fm, Luff > -70)))
        if np.logical_or(np.size(np.nonzero(20*np.log10(np.abs(Hf[idx])) > Luff[idx])) != 0,
                         np.size(np.nonzero(20*np.log10(np.abs(Hf[idx])) < Llff[idx])) != 0):
            print("\nWarning: filter shape does not conform with BS EN "
                  "61620:2014 Class 1 accuracy limits. Resample data to a "
                  "different sampling frequency or select a different filter "
                  "order.")
            plt.figure()
            plt.plot(ff, 20*np.log10(abs(Hf)), label="Filter")
            plt.plot(ff, Luff, label="IEC Upper limit")
            plt.plot(ff, Llff, label="IEC lower limit")
            G10 = 10**(3/10)
            octRatio = G10**(0.5/n)
            plt.xlim((fm/octRatio/1.5, fm*octRatio*1.5))
            plt.ylim((-100, 10))
            plt.xlabel("Frequency, Hz")
            plt.ylabel("Filter response, dB")
            plt.legend(loc=4)

    # TODO
    # # if necessary, downsample signal for filter
    # if fs_filt != fs:
    #     x_filt = resample_poly(x, up=1, down=int(fs/fs_filt),
    #                            padtype='line', axis=axis)
    # else:
    #     x_filt = x

    if fwd_bwd:
        # filter signal using forward-backward processing
        # generate fbw pad length based on Boore, 2005
        padN = np.min([x.shape[axis]-1, int(x.shape[axis]*1.5*order/2/f2.item())])
        y_filt = sosfiltfilt(sos, x, axis=axis, padtype='constant',
                             padlen=padN)
    else:
        y_filt = sosfilt(sos, x, axis=axis)  # filter signal

    # # resample to original frequency
    # if fs_filt != fs:
    #     y = resample_poly(y_filt, up=int(fs/fs_filt), down=1,
    #                       padtype='line', axis=axis)
    # else:
    y = y_filt

    return y


def noct2spec(Xnoct, fmoct, xintype='dB', spectype='rms', axis=0):
    """
    Converts input fractional octave band spectrum to an interpolated
    narrowband equivalent.

    Parameters
    ----------
    Xnoct : 1D or 2D array
            Input fractional 1/N octave band spectral values, either in dB or
            linear units.
    fmoct : 1D array
            Input band mid frequencies for xnoct.
            Must have the same number of frequencies as xnoct.
    Xintype: string, 'dB' or 'lin', optional (default: 'dB')
             The form of Xnoct, either dB ('dB') or linear units ('lin').
    spectype : string, 'rms', 'psd' or 'frf', optional (default: 'rms')
               The type of output spectrum, either RMS amplitude ('rms'), power
               spectral density ('psd'), or frequency response function ('frf').
    axis : integer, 0 or 1, optional (default: 0)
           Axis along which the input spectra are defined (0=rows, 1=columns).

    Returns
    -------
    Xspec : 1D or 2D array
        Output narrowband equivalent spectrum.
    fspec : TYPE
        DESCRIPTION.

    """

    #TODO return Xspec, fspec
