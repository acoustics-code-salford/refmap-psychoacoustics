# -*- coding: utf-8 -*-
# %% Preamble
"""
acousticSHMVerification.py
--------------

Acoustic signal analysis subfunctions for implementing the Sottek Hearing
Model, as defined in the ECMA-418-2 standard (currently 2025).

Requirements
------------
soundfile
scipy
refmap-psychoacoustics (metrics.ecma418_2)

Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford

Date created: 22/07/2025
Date last modified: 22/07/2025
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
# %% Import block
import os
import soundfile as sf
from scipy.io import savemat
from src.py.metrics.ecma418_2 import (acousticSHMTonality, acousticSHMLoudness,
                                      acousticSHMLoudnessFromComponent,
                                      acousticSHMRoughness)


# %% acousticSHMVerification
def acousticSHMVerification():
    # %% Load paths (assumes root directory is refmap-psychoacoustics)
    rootPath = os.getcwd()
    inputPath = os.path.join(rootPath, 'validation', 'ECMA-418-2', 'audio')
    outputPath = os.path.join(rootPath, 'validation', 'ECMA-418-2', 'py')

    # %% Import audio

    p, sampleRateIn = sf.read(os.path.join(inputPath, 'BusyStreet1_0530-0600.wav'))

    # %% Calculate metric values

    tonalitySHM = acousticSHMTonality(p, sampleRateIn, axisN=0,
                                      soundField='freeFrontal', waitBar=True,
                                      outPlot=False)

    loudnessSHM = acousticSHMLoudnessFromComponent(tonalitySHM['specTonalLoudness'],
                                                   tonalitySHM['specNoiseLoudness'],
                                                   outPlot=False, binaural=True)

    roughnessSHM = acousticSHMRoughness(p, sampleRateIn, axisN=0,
                                        soundField='freeFrontal', waitBar=True,
                                        outPlot=False, binaural=True)

    # %% store outputs as .mat files
    savemat(os.path.join(outputPath, "BusyStreet1_0530-0600_tonalitySHM_Py.mat"),
            tonalitySHM, appendmat=False, do_compression=True,
            oned_as='column')

    savemat(os.path.join(outputPath, "BusyStreet1_0530-0600_loudnessSHM_Py.mat"),
            loudnessSHM, appendmat=False, do_compression=True,
            oned_as='column')

    savemat(os.path.join(outputPath, "BusyStreet1_0530-0600_roughnessSHM_Py.mat"),
            roughnessSHM, appendmat=False, do_compression=True,
            oned_as='column')

