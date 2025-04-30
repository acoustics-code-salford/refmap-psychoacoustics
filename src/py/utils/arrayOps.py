# -*- coding: utf-8 -*-
# %% Preamble
"""
arrays.py
--------------

Array manipulations and operations

Requirements
------------
numpy

Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford

Date created: 30/04/2025
Date last modified: 30/04/2025
Python version: 3.11

Copyright statement: This file and code is part of work undertaken within
the RefMap project (www.refmap.eu), and is subject to licence as detailed
in the code repository
(https://github.com/acoustics-code-salford/refmap-psychoacoustics)

As per the licensing information, please be aware that this code is WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

Checked by:
Date last checked:

"""

# %% Imports and parameter setup
import numpy as np


# %% dimensional
def dimensional(ndArray, targetDim=2, where='first'):
    """
    Return array increased by dimensions depending on difference between
    targetDim and len(ndarray.shape), placed either 'first' or 'last' according
    to where keyword argument

    Parameters
    ----------
    ndArray : array_like
        Input array.
    targetDim : TYPE, optional
        The target number of dimensions. The default is 2.
    where : keyword string or corresponding integer (0, -1), optional
        Where the added dimensions are to be placed.
        The default is 'first' (0). The altgernative is 'last' (-1).

    Returns
    -------
    targArray : nD array
                Output numpy array increased by dimensions.

    """

    # check input type and try to convert to numpy array
    if type(ndArray) is np.ndarray:
        pass
    else:
        try:
            ndArray = np.array(ndArray)
        except TypeError:
            raise TypeError("Input ndArray must be an array-like object.")

    # assign dimensions to add
    dimsToAdd = targetDim - ndArray.ndim

    # add number of dimensions
    targArray = ndArray
    if dimsToAdd <= 0:
        return targArray
    else:
        if where == 'first' or where == 0:
            for ii in range(dimsToAdd):
                targArray = np.expand_dims(targArray, 0)
        elif where == 'last' or where == -1:
            for ii in range(dimsToAdd):
                targArray = np.expand_dims(targArray, -1)
        else:
            raise ValueError("Input argument 'where' must either have string values 'first' or 'last', or integer values 0 or -1")

    return targArray
