# -*- coding: utf-8 -*-
"""
fft_tools.py
------------

Module provides tools commonly used in FFT spectral analysis

Requirements
------------
numpy (1.23.4)
scipy (1.9.3)

Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford
 
Date created: 26/10/2023
Date last modified: 14/07/2024
Python version: 3.11

Copyright statement: This file and code is part of work undertaken within
the RefMap project (www.refmap.eu), and is subject to licence as detailed
in the code repository
(https://github.com/acoustics-code-salford/refmap-psychoacoustics)

Checked by:
Date last checked:

"""

from scipy.fft import (fft, fftfreq, fftshift, rfft, rfftfreq)
from scipy.signal import (get_window)
import numpy as np


def fft_spec_Welch(x, fs, nfft, over=0.5, axis=0, win_type='hann', twosided=False):
    '''
    Return the one- or two-sided FFT magnitude spectrum of x with associated
    frequencies, using Welch's method of averaged periodograms - outputs can
    include any of the following spectrum types:
        Linear peak amplitude
        Linear RMS amplitude
        Power (auto) spectral density
        Autospectrum (mean-square amplitude, ie squared linear RMS amplitude)

    One-sided transform output will be nfft/2 in length if nfft is even,
    or nfft/2 + 1 if nfft is odd. Two-sided spectrum will be nfft in length.

    Inputs
    ------
    x : 1D or 2D array of vectors
        contains the input signals (in rows or columns specified by 'axis',
        if 2D)
    fs : integer
        the signal sampling frequency
    nfft : integer
        the length of the transform block, related to the frequency interval
        'df' by nfft = fs/df
    over : float between 0.0 and (<)1.0
        the overlap to use for windowed averaging (only applies if signal
        length is at least two overlapped window lengths)
    axis : integer
        the axis of the input array along which to apply the transform
    win_type : string or None
        indicates the type of window to apply to the block
    twosided : Boolean
               indicates the type of spectrum to output (one- or two-sided)
   
    Returns
    -------
    fftf : 1D array
        contains the FFT spectral line frequencies, Hz
    Ben : float
        the normalised effective (equivalent) noise bandwidth for the window function applied
    Lx: 1D or 2D array
        contains the linear peak magnitude spectrum of x (units of amplitude)
    LxRMS : 1D or 2D array
        contains the linear RMS magnitude spectrum of x  (units of amplitude)
    Gxx : 1D or 2D array
        contains the power spectral density of x (units of amplitude**2/Hz)
    Axx : 1D or 2D array
        contains the autospectrum of x (ie the square of LxRMS, units of
        amplitude**2)

    Assumptions
    -----------
    x is a 1D signal or 2D array of signals
    '''
    if nfft > x.shape[axis]:
        raise ValueError("The FFT length must be no longer than the input signal.")
    
    if over < 0 or over >= 1:
        raise ValueError("The overlap proportion must have a non-negative value, 0 <= overlap <1.")
    

    if win_type is None:
        win_type = 'boxcar'

    dt = 1/fs  # timestep, s
    df = fs/nfft  # frequency interval, Hz

    window = get_window(win_type, nfft)  # generate window of fft block length

    Awin = nfft/np.sum(window)  # amplitude correction for windowing
    Be = df*nfft*np.sum(window**2)/np.sum(window)**2  # window equivalent noise bandwidth
    Ben = Be/df  # window normalised equivalent noise bandwidth

    if np.mod(nfft, 2) == 0:
        odd = 0
    else:
        odd = 1

    if twosided:
        # two-sided range of positive frequencies for spectral lines
        fftf = fftshift(fftfreq(nfft, dt))
    else:
        # one-sided range of positive frequencies for spectral lines
        fftf = fftfreq(nfft, dt)[:int(nfft//2) + odd]

    iiEnd = fftf.shape[0]  # index for the spectral size

    # transpose array of row vectors to array of column vectors
    if axis == 1:
        x = np.transpose(x)

    axis2 = 0  # dummy axis variable

    if len(x.shape) == 1:
        Axx = np.zeros(fft(x[0:nfft]).shape)
    else:
        Axx = np.zeros(fft(x[0:nfft, :], axis=0).shape)

    block = 0  # initialise block counter
    blstart = 0  # block start index
    blend = int(nfft)  # block end index

    # loop over blocks and accumulate unscaled spectra
    while blend <= x.shape[axis2]:

        if len(x.shape) == 1:
            xt = window*x[blstart:blend]
            XX = np.real(np.conj(fft(xt))*fft(xt))
        else:
            win = np.moveaxis(np.broadcast_to(window,
                                              (x.shape[1],
                                               np.size(window))), -1, 0)
            xt = win*x[blstart:blend, :]
            XX = np.real(np.conj(fft(xt, axis=0))*fft(xt, axis=0))

        Axx += XX
        block += 1
        blstart = int(block*nfft*(1 - over))
        blend = int(blstart + nfft)

    Axx = Axx/block  # average the autospectrum over blocks
    Axx = (Awin/nfft)**2*Axx  # window amplitude and FFT correction

    if twosided is False:
        Axx = Axx[:iiEnd]
        if len(x.shape) == 1:
            Axx[1:] = 2*Axx[1:]  # compensate for single-sided autospectrum
        else:
            Axx[1:, :] = 2*Axx[1:, :]  # compensate for single-sided autospectrum
    else:
        Axx = fftshift(Axx)

    # transpose output array back to input row vector form
    if axis == 1:
        Axx = np.transpose(Axx)

    LxRMS = np.sqrt(Axx)  # linear RMS spectrum
    Lx = LxRMS*np.sqrt(2)  # linear peak spectrum
    Gxx = Axx/Be  # corrects autospectrum for effective noise bandwidth to give PSD

    if __name__ == '__main__' and block > 1:
        print("\nFFT processing completed using " + str(int(block)) +
              " windowed averages and " + str(int(over*100)) +
              "% window overlap")

    return fftf, Ben, Lx, LxRMS, Gxx, Axx


def rfft_spec_Welch(x, fs, nfft, over=0.5, axis=0, win_type='hann'):
    '''
    Return the one-sided FFT magnitude spectrum of a real-valued signal x,
    with associated frequencies, using Welch's method of averaged periodograms
    - outputs can include any or all of the following spectrum types:
        Linear peak amplitude
        Linear RMS amplitude
        Power (auto) spectral density
        Autospectrum (mean-square amplitude, ie squared linear RMS amplitude)

    Transform output will be nfft//2 + 1 in length

    Inputs
    ------
    x : 1D or 2D array of vectors
        contains the input signals (in rows or columns specified by 'axis',
        if 2D)
    fs : integer
        the signal sampling frequency
    nfft : integer
        the length of the transform block, related to the frequency interval
        'df' by nfft = fs/df
    over : float between 0.0 and (<)1.0
        the overlap to use for windowed averaging (only applies if signal
        length is at least two overlapped window lengths)
    axis : integer
        the axis of the input array along which to apply the transform
    win_type : string or None
        indicates the type of window to apply to the block
   
    Returns
    -------
    fftf : 1D array
        contains the FFT spectral line frequencies, Hz
    Ben : float
        the normalised effective (equivalent) noise bandwidth for the window function applied
    Lx: 1D or 2D array
        contains the linear peak magnitude spectrum of x (units of amplitude)
    LxRMS : 1D or 2D array
        contains the linear RMS magnitude spectrum of x  (units of amplitude)
    Gxx : 1D or 2D array
        contains the power spectral density of x (units of amplitude**2/Hz)
    Axx : 1D or 2D array
        contains the autospectrum of x (ie the square of LxRMS, units of
        amplitude**2)

    Assumptions
    -----------
    x is a real-valued 1D signal or 2D array of signals
    '''
    if nfft > x.shape[axis]:
        raise ValueError("The FFT length must be no longer than the input signal.")
    
    if over < 0 or over >= 1:
        raise ValueError("The overlap proportion must have a non-negative value, 0 <= overlap <1.")
    
    if win_type is None:
        win_type = 'boxcar'

    dt = 1/fs  # timestep, s
    df = fs/nfft  # frequency interval, Hz

    window = get_window(win_type, nfft)  # generate window of fft block length

    Awin = nfft/np.sum(window)  # amplitude correction for windowing
    Be = df*nfft*np.sum(window**2)/np.sum(window)**2  # window equivalent noise bandwidth
    Ben = Be/df  # window normalised equivalent noise bandwidth

    fftf = rfftfreq(nfft, dt) # one-sided range of positive frequencies for spectral lines

    # transpose array of row vectors to array of column vectors
    if axis == 1:
        x = np.transpose(x)
        
    axis2 = 0  # dummy axis variable

    if len(x.shape) == 1:
        Axx = np.zeros(rfft(x[0:nfft]).shape)
    else:
        Axx = np.zeros(rfft(x[0:nfft, :], axis=0).shape)

    block = 0  # initialise loop counter
    blstart = 0  # block start index
    blend = int(nfft)  # block end index

    while blend <= x.shape[axis2]:

        if len(x.shape) == 1:
            xt = window*x[blstart:blend]
            XX = np.real(np.conj(rfft(xt))*rfft(xt))
        else:
            win = np.moveaxis(np.broadcast_to(window, (x.shape[1], np.size(window))), -1, 0)
            xt = win*x[blstart:blend, :]
            XX = np.real(np.conj(rfft(xt, axis=0))*rfft(xt, axis=0))

        Axx += XX
        block += 1
        blstart = int(block*nfft*(1 - over))
        blend = int(blstart + nfft)

    Axx = Axx/block  # average the autospectrum over blocks
    Axx = (Awin/nfft)**2*Axx  # window amplitude and FFT correction
    if len(x.shape) == 1:
        Axx[1:] = 2*Axx[1:]  # compensate for single-sided autospectrum
    else:
        Axx[1:, :] = 2*Axx[1:, :]  # compensate for single-sided autospectrum

    # transpose output array back to input row vector form
    if axis == 1:
        Axx = np.transpose(Axx)

    LxRMS = np.sqrt(Axx)  # linear RMS spectrum
    Lx = LxRMS*np.sqrt(2)  # linear peak spectrum
    Gxx = Axx/Be  # corrects autospectrum for effective noise bandwidth to give PSD

    if __name__ == '__main__' and block > 1:
        print("\nFFT processing completed using " + str(int(block)) +
              " windowed averages and " + str(int(over*100)) + "% window overlap")

    return fftf, Ben, Lx, LxRMS, Gxx, Axx
