# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 16:21:31 2023

@author: m_lot
"""

import numpy as np
import pandas as pd
import mosqito
from scipy import optimize


def sqm_tvar(x_Pa, fs, pcex=5):
    '''
    Return time-varying sound quality metric values for an audio input file.

    '''


# Bolt, Beranek & Newman detectability model
# ------------------------------------------

def detectBBN():
    """
    Returns
    -------
    None.

    """


def _detectBBNEfficiency():
    """


    Returns
    -------
    TYPE
        DESCRIPTION.

    """
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