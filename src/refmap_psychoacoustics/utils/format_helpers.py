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
    return round(val + 10**(-len(str(val)) - 1), digits)


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
    
    crit = 1/10**digits

    if abs(val) < crit and floor:
        if val > 0:
            val_string = "<" + str(copysign(crit, val))
        else:
            val_string = ">" + str(copysign(crit, val))

    else:
        valrnd = round_trad(val, digits)
        val_string = f'{valrnd:f}'.rstrip("0")
        dec = int(val_string.split('.')[-1])
        zs = len(val_string.split('.')[-1]) - len(val_string.split('.')[-1].lstrip("0"))
        nzeros = digits - len(str(dec)) - zs
        if nzeros > 0:
            val_string = val_string + nzeros*"0"

    return(val_string)