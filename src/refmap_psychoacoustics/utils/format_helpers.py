# -*- coding: utf-8 -*-
"""
format_helpers.py
-----------------

Formatting utility functions for use in the RefMap project.

Functions
---------

display_round(val, digits=3, floor=True)
    Return a string representation of a value rounded to a specified number of decimal places.

round_trad(val, digits=3)
    Return a rounded value to a specified number of decimal places using the traditional approach.

Requirements
------------


Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford
 
Date created: 10/03/2025
Date last modified: 10/03/2025
Python version: 3.11

Copyright statement: This file and code is part of work undertaken within
the RefMap project (www.refmap.eu), and is subject to licence as detailed
in the code repository
(https://github.com/acoustics-code-salford/refmap-psychoacoustics)

Checked by:
Date last checked:

"""

from math import copysign
import numpy as np

def round_trad(val, digits=3):
    """
    Return a rounded value to a specified number of decimal places using the traditional approach.

    Parameters
    ----------
    
    val : float
          The value to be rounded.

    digits : int
             The number of decimal places to round to (default=3).
    
    Returns
    -------

     : float
          The rounded value.

    """
    # class check for numpy arrays
    if hasattr(val, 'dtype') or isinstance(val, (list, tuple)):
        arr = np.asarray(val, dtype=float)
        # compute a per-element epsilon
        eps = np.vectorize(lambda v: 10**(-len(str(v)) - 1))(arr)
        out = np.round(arr + eps, digits)
    else:
        out = round(val + 10**(-len(str(val)) - 1), digits)

    return out


def display_round(val, digits=3, floor=True):
    """
    Return a string representation of a value rounded to a specified number of decimal places.

    Parameters
    ----------

    val : float
          The value to be rounded.
    
    digits : int
             The number of decimal places to round to (default=3).

    floor : bool
            If True, the string will represent the value as "<" or ">" signed 1/10**digits

    Returns
    -------

    val_string : str
                The rounded value as a string.
        
    """

    # check if array-like and apply function recursively
    if isinstance(val, (list, tuple, np.ndarray)):
        return np.array(
            [display_round(v, digits=digits, floor=floor) for v in val],
            dtype=object,
        )
    
    crit = 1/10**digits

    if abs(val) < crit and floor:
        if val > 0:
            val_string = "<" + str(copysign(crit, val))
        else:
            val_string = str(copysign(crit, val))

    else:
        valrnd = round_trad(val, digits)
        val_string = f'{valrnd:.{digits}f}'

    return(val_string)