# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 16:21:31 2023

@author: m_lot
"""

import numpy as np
import pandas as pd
from scipy import (optimize, signal)
from src.py.dsp import noct
from math import gcd
import bottleneck as bn


# Bolt, Beranek & Newman detectability model as developed further by NASA
# -----------------------------------------------------------------------------

def detectBBNNASA(signalTarget, sampleRateTarget, axisTarget,
                  signalMasker, sampleRateMasker, axisMasker):
    """
    Returns
    -------
    
    
    Assumptions
    -----------
    Both input signals are calibrated to absolute acoustic pressure (Pascals),
    or otherwise, equal scaling.

    """
    
    if axisTarget == 1:
        signalTarget = np.transpose(signalTarget)
    
    if axisMasker == 1:
        signalMasker = np.transpose(signalMasker)
        
    # check input signal sampling frequencies are equal, otherwise resample to
    # match higher rate
    if sampleRateMasker > sampleRateTarget:
        # upsampled sampling frequency
        sampleRate = sampleRateMasker
        up = int(sampleRate/gcd(sampleRate, sampleRateTarget))
        down = int(sampleRateTarget/gcd(sampleRate, sampleRateTarget))
        signalTarget = signal.resample_poly(signalTarget, up, down, padtype='line')
        
    elif sampleRateTarget > sampleRateMasker:
        sampleRate = sampleRateTarget
        up = int(sampleRate/gcd(sampleRate, sampleRateMasker))
        down = int(sampleRateMasker/gcd(sampleRate, sampleRateMasker))
        signalMasker = signal.resample_poly(signalMasker, up, down, padtype='line')
    else:
        sampleRate = sampleRateTarget
        
    timeStep = 0.5
    timeSteps = int(0.5*sampleRate)
    
    # check signal lengths match
    if signalTarget.shape[0] != signalMasker.shape[0]:
        raise ValueError("The lengths of the input signals do not match.")
    
    # check number of channels match
    if signalTarget.shape[1] != signalMasker.shape[1]:
        raise ValueError("The numbers of channels in the input signals do not match.")
    else:
        numChans = signalTarget.shape[1]

    # it should be quicker to use a STFT to obtain spectra that can then be
    # integrated over frequency, however, attempts to use SciPy
    # signal.ShortTimeFFT with 'psd' scaling did not give outputs scaled in a
    # manner commensurate with expected results. Hence, time domain filters are
    # applied here in a loop, to ensure consistency with the NASA method,
    # presumed to be operated using MATLAB poctave

    fm, f1, f2 = noct.noctf(20, min(sampleRate/2.4, 20000), 3)
    # TODO finish algorithm
    

def _detectBBNEfficiency():
    """


    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    # curve from Figure 5 of Fidell & Horonjeff, 1982, A graphic method for
    # predicting audibility of noise sources
    detectEfficiencyData = pd.DataFrame(data=np.array([[32.3, 39.8, 51.2, 66.2,
                                                       89.1, 119.9, 161.4,
                                                       217.2, 292.3, 393.4,
                                                       529.5, 654.1, 800.5,
                                                       979.9, 1290.8, 1737.3,
                                                       2338.2, 3146.9, 4235.4,
                                                       5700.3, 7671.9, 10325.5,
                                                       13441.1],
                                                      [0.134, 0.163, 0.203,
                                                       0.247, 0.293, 0.334,
                                                       0.369, 0.394, 0.411,
                                                       0.422, 0.429, 0.430,
                                                       0.431, 0.430, 0.424,
                                                       0.414, 0.397, 0.373,
                                                       0.343, 0.310, 0.276,
                                                       0.243, 0.213]]).T,
                                        columns=["Hz", "eta"])

    p0 = [19, 0.02, 1, 0.14, 140, -0.9]

    detectEffCurveFit, _ = optimize.curve_fit(f=_detectEfficiencyFunc,
                                              xdata=detectEfficiencyData['Hz'],
                                              ydata=detectEfficiencyData['eta'],
                                              p0=p0, maxfev=1000000,
                                              bounds=(-10, 10000))
    detectEffCurveTestFit = _detectEfficiencyFunc(detectEfficiencyData['Hz'],
                                                  detectEffCurveFit[0],
                                                  detectEffCurveFit[1],
                                                  detectEffCurveFit[2],
                                                  detectEffCurveFit[3],
                                                  detectEffCurveFit[4],
                                                  detectEffCurveFit[5])
    
    return detectEffCurveFit

def _detectEfficiencyFunc(freq, a, b, c, d, e, f):
    return (freq/800)**a/(b + c*(freq/800)**d)**e + f