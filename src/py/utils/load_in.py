# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 18:11:43 2023

@author: m_lot
"""

import os
import numpy as np
import librosa
from src.py.dsp.filterFuncs import A_weight_T


def audio_caltnc(filepath, cal_val, tnc=[0, 0], cal_tnc='pre', cal_type=None):
    '''
    Returns calibrated, truncated (optional) signal from audio input file.
    Inputs
    ------
    filepath : string
               the path to the input audio file
    tnc : vector or 2D array of floats, length 2
          the truncation values (in seconds) to be applied to the start or
          end of the input file - the first value in tnc applies to the start
          and the second value to the end
    cal_val : double
              the value of the calibration to be applied to the file
    cal_tnc : keyword string (default: 'pre')
              specifies whether calibration of the input file should take
              place pre- or post-truncation operation
              'pre':apply calibration before file truncation
              'post': apply calibration to the truncated file
    cal_type : keyword string (default: None)
               specifies the type of calibration for the input file (ie, what
               property cal_val represents)
    
    Returns
    -------
    x_Pa : vector or 2D array
            the audio signal sample data calibrated to units of Pascals
       
    fs : positive integer
         the sample frequency of the audio signal
       
    Ownership and Quality Assurance
       
    Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
    Institution: University of Salford
     
    Date created: 21/07/2023
    Date last modified: 30/04/2025
    Python version: 3.10.9
       
    Copyright statement: This file and code is part of work undertaken within
    the RefMap project (www.refmap.eu), and is subject to licence as detailed
    in the code repository
    (https://github.com/acoustics-code-salford/refmap-psychoacoustics)
       
    Checked by:
    Date last checked:
       
    '''
    ## Input, truncation and calibration

    xt, fs = librosa.load(filepath, sr=None, mono=False)
    xt = np.transpose(xt)

    # signal calibration block (if pre-truncation)
    if cal_tnc == 'pre':
        match cal_type.upper():
            case 'LAEQ':
                xA_temp = A_weight_T(xt, fs)
                xA_rms = np.sqrt(np.mean(xA_temp**2, axis=0))
                cal_valPa = 2e-5*10**(cal_val/20)
                xn = xt*(cal_valPa/xA_rms)
                del xA_temp
            case 'LAE':
                xA_temp = A_weight_T(xt, fs)
                xA_rss = np.sqrt(xA_temp.shape[0]/fs)*np.sqrt(np.mean(xA_temp**2, axis=0))
                cal_valPa = 2e-5*10**(cal_val/20)
                xn = xt*(cal_valPa/xA_rss)
                del xA_temp
            case 'LPEAK':
                xpeak = np.max(np.abs(xt), axis=0)
                cal_valPa = 2e-5*10**(cal_val/20)
                xn = xt*(cal_valPa/xpeak)
    else:
        xn = xt


    # signal truncation block
    # start truncation
    if tnc[0] != 0:
        if tnc[0]*fs < np.shape(xn)[0]:
            tnc_s = int(np.floor(tnc[0]*fs))
            if len(np.shape(xn)) > 1:
                xn = xn[tnc_s:, :]
            else:
                xn = xn[tnc_s:]
        else:
            raise ValueError(f"Truncation value {tnc[0]} is larger than signal length {np.shape(xn, 0)}")

    # end truncation
    if tnc[1] != 0:
        if tnc[1]*fs < np.shape(xn)[0]:
            tnc_e = int(np.floor(tnc[1]*fs))
            if len(np.shape(xn)) > 1:
                xn = xn[:-tnc_e, :]
            else:
                xn = xn[:-tnc_e]
        else:
            raise ValueError(f"Truncation value {tnc[1]} is larger than signal length {np.shape(xn, 0)}")
        


    # signal calibration block (if post-truncation)
    if cal_tnc == 'post':
        match cal_type.upper(): # calibrate to type
            case 'LAEQ':
                xA_temp = A_weight_T(xn, fs)
                xA_rms = np.sqrt(np.mean(xA_temp**2, axis=0))
                cal_valPa = 2e-5*10**(cal_val/20)
                x_Pa = xn*(cal_valPa/xA_rms)
                del xA_temp
            case 'LAE':
                xA_temp = A_weight_T(xn, fs)
                xA_rss = np.sqrt(xA_temp.shape[0]/fs)*np.sqrt(np.mean(xA_temp**2, axis=0))
                cal_valPa = 2e-5*10**(cal_val/20)
                x_Pa = xn*(cal_valPa/xA_rss)
                del xA_temp
            case 'LPEAK':
                xpeak = np.max(np.abs(xn), axis=0)
                cal_valPa = 2e-5*10**(cal_val/20)
                x_Pa = xn*(cal_valPa/xpeak)
        
    else:
        x_Pa = xn

    
    return x_Pa, fs
