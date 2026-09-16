# -*- coding: utf-8 -*-
# %% Preamble
"""
data_helpers.py
--------------

Data transformation and manipulation functions

Functions
---------

ray_polygon_intersection(direction, polygon)
    Calculate the intersection of a ray with a convex polygon enveloping the origin.

iso_transform_circ(x, y)
    Transform input coordinate points from a locus bounded by a circle to the 
    domain bounded by an octagon as defined by the ISO 12913-3:2025 Annex A
    transform functions.

Requirements
------------
numpy
matplotlib

Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford

Date created: 22/04/2026
Date last modified: 22/04/2026
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
from matplotlib import pyplot as plt

# %% Function definitions

# Calculate the intersection of a ray with a convex polygon enveloping the origin
def ray_polygon_intersection(direction, polygon):
    """
    Return a scaling factor t such that the ray defined by the
    input direction intersects with an input convex polygon
    enveloping the origin.

    """

    dx, dy = direction
    t_min = np.inf

    for i in range(len(polygon)):
        p1 = polygon[i]
        p2 = polygon[(i + 1) % len(polygon)]
        
        edge = p2 - p1
        A = np.array([[dx, -edge[0]],
                      [dy, -edge[1]]])
        
        if abs(np.linalg.det(A)) < 1e-12:
            continue
        
        t, u = np.linalg.solve(A, p1)
        
        if t > 0 and 0 <= u <= 1:
            t_min = min(t_min, t)
    
    return t_min


def iso_transform_circ(x, y):
    """
    Return transformed coordinates (x', y') for input coordinates (x, y) by calculating
    the polar ray intersection with the corresponding octagon boundary as defined by the
    ISO 12913-3:2025 Annex A transform functions.

    Inputs
    ------
    x: float
        X-coordinate of the input point.
    
    y: float
        Y-coordinate of the input point.

    Outputs
    -------
    x': float
        X-coordinate of the transformed point.
       
    y': float
        Y-coordinate of the transformed point.

    """

    r = np.hypot(x, y)
    
    if r == 0:
        return 0.0, 0.0

    direction = np.array([x, y])
    direction = direction / r
    
    # Define the ISO octagon vertices
    a = 1.0
    b = np.sqrt(2) - 1.0

    octagon = np.array([
        [-a, -b],
        [-b, -a],
        [ b, -a],
        [ a, -b],
        [ a,  b],
        [ b,  a],
        [-b,  a],
        [-a,  b]
    ])

    t_boundary = ray_polygon_intersection(direction, octagon)
    
    return t_boundary*x, t_boundary*y