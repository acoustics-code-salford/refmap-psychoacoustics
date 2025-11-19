"""
psych_annoy.py
----------

Psychoacoustic annoyance metric implementations:

widmannPA : Widmann (1992)
morePA : More (2010)
willemsenPA : Willemsen et al. (2010)
diPA : Di et al. (2016)
torijaPA : Torija et al. (2022)

Requirements
------------
numpy

References
----------
See individual function docstrings for references.

Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford

Date created: 15/11/2025
Date last modified: 19/11/2025
Python version: 3.11

Copyright statement: This file and code is part of work undertaken within
the RefMap project (www.refmap.eu), and is subject to licence as detailed
in the code repository
(https://github.com/acoustics-code-salford/refmap-psychoacoustics)

As per the licensing information, please be aware that this code is WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

"""

# %% Import block
import numpy as np


# %% widmannPA
def widmannPA(loud5ExZwicker, sharp5ExWidmann, rough5ExFastl, fluct5ExFastl):
    """
    Returns Widmann (1992) psychoacoustic annoyance metric based on time-aggregated sound quality
    metrics loudness, sharpness, roughness and fluctuation strength. The sound quality metrics are
    assumed to be 5%-exceeded (95th percentile) values.

    Parameters
    ----------
    loud5ExZwicker : float or array
        Zwicker loudness values as defined in ISO 532-1:2017, 5%-exceeded (95th percentile), sone.

    sharp5ExWidmann : float or array
        Widmann-weighted sharpness values as defined in DIN 45692:2009 and Widmann (1992)
        (using Zwicker loudness input), 5%-exceeded (95th percentile), acum.

    rough5ExFastl : float or array
        Fastl roughness values as defined in Fastl (1977) and Fastl & Zwicker (2007),
        5%-exceeded (95th percentile), asper.

    fluct5ExFastl : float or array
        Fastl fluctuation strength values as defined in Fastl (1982) and Fastl & Zwicker (2007),
        5%-exceeded (95th percentile), vacil.

    Returns
    -------
    psych_annoy : float or array
        Widmann (1992) psychoacoustic annoyance metric, au.

    Assumptions
    -----------
    All input metrics are same shape (scalar or array), non-negative, and represent time-aggregated
    5%-exceeded (95th percentile) values.

    References
    ----------
    ISO 532-1:2017. Acoustics — Method for calculating loudness — Part 1: Zwicker method.

    DIN 45692:2009. Messtechnische Simulation der Hörempfindung Schärfe (Measurement technique for
    the simulation of the auditory sensation of sharpness).

    Fastl, H. (1977). Roughness and temporal masking patterns of sinusoidally amplitude
    modulated broadband noise. Psychophysics and Physiology of Hearing, Keele, UK , 12–16 April 1977.

    Fastl, H. (1982). Fluctuation strength and temporal masking patterns of amplitude-modulated
    broadband noise. Hearing Research, 8(1), 59–69.

    Fastl, H. & Zwicker, E. (2007). Psychoacoustics: Facts and models (3rd ed.). Springer.

    Widmann, U. (1992). Ein Modell der Psychoakustischen Lästigkeit von Schallen und seine Anwendung
    in der Praxis der Lärmbeurteilung (A model of the psychoacoustic annoyance of sounds and its
    application in noise assessment practice) [Doctoral thesis, Technische Universität München
    (Technical University of Munich)]

    """
    wS = (sharp5ExWidmann - 1.75)*0.25*np.log10(loud5ExZwicker + 10)
    wS[wS <= 1.75] = 0
    wMod = 2.18/(loud5ExZwicker)**0.4*(0.4*fluct5ExFastl + 0.6*rough5ExFastl)

    psych_annoy = loud5ExZwicker*(1 + np.sqrt(wS**2 + wMod**2))

    return psych_annoy


# %% morePA
def morePA(loud5ExZwicker, sharp5ExMore, rough5ExFastl, fluct5ExFastl, tonal5ExAures):
    """
    Returns More et al. (2010) psychoacoustic annoyance metric based on time-aggregated sound quality
    metrics loudness, sharpness, roughness, fluctuation strength and tonality. The sound quality
    metrics are assumed to be 5%-exceeded (95th percentile) values.

    Parameters
    ----------
    loud5ExZwicker : float or array
        Zwicker loudness values as defined in ISO 532-1:2017, 5%-exceeded (95th percentile), sone.

    sharp5ExMore : float or array
        More-weighted sharpness values as defined in More (2011, using Zwicker loudness input),
        5%-exceeded (95th percentile), acum. Note that More's sharpness weighting differs slightly from
        the von Bismarck (1974) weighting that formed the weighting republished in Zwicker & Fastl (1999),
        which More attempted to model using an exponential function. In lieu of More-weighted sharpness values,
        von Bismarck-weighted sharpness values should be used, as the nearest approximation.

    rough5ExFastl : float or array
        Fastl roughness values as defined in Fastl (1977) and Zwicker & Fastl (1999),
        5%-exceeded (95th percentile), asper.

    fluct5ExFastl : float or array
        Fastl fluctuation strength values as defined in Fastl (1982) and Zwicker & Fastl (1999),
        5%-exceeded (95th percentile), vacil.
    
    tonal5ExAures : float or array
        Aures tonality values as defined in Aures (1985), 5%-exceeded (95th percentile), tu.

    Returns
    -------
    psych_annoy : float or array
        More (2010) psychoacoustic annoyance metric, au.

    Assumptions
    -----------
    All input metrics are same shape (scalar or array), non-negative, and represent time-aggregated
    5%-exceeded (95th percentile) values.

    References
    ----------
    ISO 532-1:2017. Acoustics — Method for calculating loudness — Part 1: Zwicker method.

    Fastl, H. (1977). Roughness and temporal masking patterns of sinusoidally amplitude
    modulated broadband noise. Psychophysics and Physiology of Hearing, Keele, UK , 12–16 April 1977.

    Fastl, H. (1982). Fluctuation strength and temporal masking patterns of amplitude-modulated
    broadband noise. Hearing Research, 8(1), 59–69.

    More, S. & Davies, P. (2010). Development of a model to predict annoyance caused by aircraft noise.
    In: Proceedings of Inter-noise 2010, Lisbon, Portugal, 13-16 June 2010.

    More, S. R. (2011). Aircraft noise characteristics and metrics. Doctoral thesis, Purdue University.

    Aures, W. (1985). Berechnungsverfahren für den sensorischen Wohlklang beliebiger Schallsignale
    (Calculation method for the sensory euphony of arbitrary sound signals). Acustica, 59(2), 130–141.

    von Bismarck, G. (1974). Sharpness as an attribute of the timbre of steady sounds.
    Acustica, 30(3), 159–172.

    Zwicker, E. & Fastl, H. (1999). Psychoacoustics: Facts and models (2nd ed.). Springer.

    """
    wS = (sharp5ExMore - 1.75)*0.25*np.log10(loud5ExZwicker + 10)
    wS[wS <= 1.75] = 0
    wMod = 2.18/(loud5ExZwicker)**0.4*(0.4*fluct5ExFastl + 0.6*rough5ExFastl)
    wT = (1 - np.exp(-0.29*loud5ExZwicker))*(1 - np.exp(-5.49*tonal5ExAures))

    sqm = np.maximum(0, -0.16 + 11.48*(wS**2) + 0.84*(wMod**2) + 1.25*(wT**2))

    psych_annoy = loud5ExZwicker*(1 + np.sqrt(sqm))

    return psych_annoy


# %% willemsenPA
def willemsenPA(loud5ExZwicker, sharpMedAures, rough5ExHEAD, impulsAvgWillemsen):
    """
    Returns Willemsen et al. (2010) psychoacoustic annoyance metric based on time-aggregated sound quality
    metrics loudness, sharpness, roughness, fluctuation strength, tonality and impulsive loudness. The sound
    quality metrics have various time aggregations determined from experimental data.

    Parameters
    ----------
    loud5ExZwicker : float or array
        Zwicker loudness values as defined in ISO 532-1:2017, 5%-exceeded (95th percentile), sone.

    sharpMedAures : float or array
        Median Aures sharpness values as defined in Aures (1985), 5%-exceeded (95th percentile), acum.

    rough5ExHEAD : float or array
        HEAD acoustics early Sottek Hearing Model roughness values as defined in ArtemiS SUITE software
        (module ACM 907), 5%-exceeded (95th percentile), asper. In lieu of HEAD roughness values, ECMA-418-2 roughness
        values may be used as the nearest approximation.

    impulsAvgWillemsen : float or array
        Willemsen (et al.) impulsive loudness values as defined in Willemsen et al. (2010),
        arithmetic mean, sone.

    Returns
    -------
    psych_annoy : float or array
        Willemsen et al. (2010) psychoacoustic annoyance metric.

    Assumptions
    -----------
    All input metrics are same shape (scalar or array), non-negative, and represent time-aggregated
    values.

    References
    ----------
    Aures, W. (1985). Berechnungsverfahren für den sensorischen Wohlklang beliebiger Schallsignale
    (Calculation method for the sensory euphony of arbitrary sound signals). Acustica, 59(2), 130–141.

    ECMA-418-2:2025. Psychoacoustic metrics for ITT equipment — Part 2 (methods for describing human
    perception based on the Sottek Hearing Model).
 
    ISO 532-1:2017. Acoustics — Method for calculating loudness — Part 1: Zwicker method.

    Willemsen, A. M. & Rao, M. D. (2010). Characterization of sound quality of impulsive sounds
    using loudness based metric. In: Proceedings of ICA 2010, Sydney, Australia, 23–27 August 2010.

    """
    psych_annoy = 0.86*loud5ExZwicker*sharpMedAures + 1.81*rough5ExHEAD + 1.24*impulsAvgWillemsen + 27.73

    return psych_annoy


# %% diPA
def diPA(loud5ExZwicker, sharp5ExAures, rough5ExHEAD, fluct5ExHEAD, tonal5ExAures):
    """
    Returns Di et al. (2016) psychoacoustic annoyance metric based on time-aggregated sound quality
    metrics loudness, sharpness, roughness, fluctuation strength and tonality. The sound quality
    metrics are assumed to be 5%-exceeded (95th percentile) values.

    Parameters
    ----------
    loud5ExZwicker : float or array
        Zwicker loudness values as defined in ISO 532-1:2017, 5%-exceeded (95th percentile), sone.

    sharp5ExAures : float or array
        Aures-weighted sharpness values as defined in Aures (1985, using Zwicker loudness input),
        5%-exceeded (95th percentile), acum.

    rough5ExHEAD : float or array
        HEAD acoustics roughness values as defined in v10 of ArtemiS SUITE software (module ACM 900),
        5%-exceeded (95th percentile), asper. In lieu of HEAD roughness values, ECMA-418-2 roughness
        values may be used as the nearest approximation.

    fluct5ExHEAD : float or array
        HEAD acoustics fluctuation strength values as defined in v10 of ArtemiS SUITE software
        (module ACM 907), 5%-exceeded (95th percentile), vacil. In lieu of HEAD fluctuation strength
        values, ECMA-418-2 fluctuation strength values may be used as the nearest approximation.
    
    tonal5ExAures : float or array
        Aures tonality values as defined in Aures (1985), 5%-exceeded (95th percentile), tu.

    Returns
    -------
    psych_annoy : float or array
        Di (2016) psychoacoustic annoyance metric, au.

    Assumptions
    -----------
    All input metrics are same shape (scalar or array), non-negative, and represent time-aggregated
    5%-exceeded (95th percentile) values.

    References
    ----------
    Di, G., Chen, X.-W., Song, K., Zhou, B. & Pei, C.-M. (2016). Improvement of Zwicker’s psychoacoustic
    annoyance model aiming at tonal noises. Applied Acoustics, 105, 164–170.

    ECMA-418-2:2025. Psychoacoustic metrics for ITT equipment — Part 2 (methods for describing human
    perception based on the Sottek Hearing Model).

    ISO 532-1:2017. Acoustics — Method for calculating loudness — Part 1: Zwicker method.

    Aures, W. (1985). Berechnungsverfahren für den sensorischen Wohlklang beliebiger Schallsignale
    (Calculation method for the sensory euphony of arbitrary sound signals). Acustica, 59(2), 130–141.

    """
    wS = (sharp5ExAures - 1.75)*0.25*np.log10(loud5ExZwicker + 10)
    wS[wS <= 1.75] = 0
    wMod = 2.18/(loud5ExZwicker)**0.4*(0.4*fluct5ExHEAD + 0.6*rough5ExHEAD)
    wT = 6.41*tonal5ExAures/loud5ExZwicker**0.52

    psych_annoy = loud5ExZwicker*(1 + np.sqrt(wS**2 + wMod**2 + wT**2))

    return psych_annoy

# %% torijaPA
def torijaPA(loud5ExZwicker, sharp5ExWidmann, rough5ExSHMOld, fluct5ExSHMOld, tonal5ExAures, impuls5ExSHMOld):
    """
    Returns Torija et al. (2022) psychoacoustic annoyance metric based on time-aggregated sound quality
    metrics loudness, sharpness, roughness, fluctuation strength, tonality and impulsiveness. The sound
    quality metrics are assumed to be 5%-exceeded (95th percentile) values.

    Parameters
    ----------
    loud5ExZwicker : float or array
        Zwicker loudness values as defined in ISO 532-1:2017, 5%-exceeded (95th percentile), sone.

    sharp5ExWidmann : float or array
        Widmann-weighted sharpness values as defined in DIN 45692:2009 and Widmann (1992)
        (using Zwicker loudness input), 5%-exceeded (95th percentile), acum.

    rough5ExSHMOld : float or array
        HEAD acoustics early Sottek Hearing Model roughness values as defined in ArtemiS SUITE software
        (module ACM 907), 5%-exceeded (95th percentile), asper. In lieu of HEAD roughness values,
        ECMA-418-2 roughness values may be used as the nearest approximation.

    fluct5ExSHMOld : float or array
        HEAD acoustics early Sottek Hearing Model fluctuation strength values as defined in ArtemiS
        SUITE software (module ACM 907), 5%-exceeded (95th percentile), vacil. In lieu of HEAD
        fluctuation strength values, ECMA-418-2 fluctuation strength values may be used as the nearest
        approximation.
    
    tonal5ExAures : float or array
        Aures tonality values as defined in Aures (1985), 5%-exceeded (95th percentile), tu.

    impuls5ExSHMOld : float or array
        HEAD acoustics early Sottek Hearing Model impulsiveness values as defined in ArtemiS SUITE
        software (module ACM 907), 5%-exceeded (95th percentile), iu. No nearest alternative metric is
        currently known.

    Returns
    -------
    psych_annoy : float or array
        More (2010) psychoacoustic annoyance metric, au.

    Assumptions
    -----------
    All input metrics are same shape (scalar or array), non-negative, and represent time-aggregated
    5%-exceeded (95th percentile) values.

    References
    ----------
    ECMA-418-2:2025. Psychoacoustic metrics for ITT equipment — Part 2 (methods for describing human
    perception based on the Sottek Hearing Model).

    ISO 532-1:2017. Acoustics — Method for calculating loudness — Part 1: Zwicker method.

    Torija, A. J., Li, Z. & Chaitanya, P. (2022). Psychoacoustic modelling of rotor noise. The Journal
    of the Acoustical Society of America, 151(3), 1804–1815.

    Widmann, U. (1992). Ein Modell der Psychoakustischen Lästigkeit von Schallen und seine Anwendung
    in der Praxis der Lärmbeurteilung (A model of the psychoacoustic annoyance of sounds and its
    application in noise assessment practice) [Doctoral thesis, Technische Universität München
    (Technical University of Munich)]

    Zwicker, E. & Fastl, H. (1999). Psychoacoustics: Facts and models (2nd ed.). Springer.

    """
    wS = (sharp5ExWidmann - 1.75)*0.25*np.log10(loud5ExZwicker + 10)
    wS[wS <= 1.75] = 0
    wMod = 2.18/(loud5ExZwicker)**0.4*(0.4*fluct5ExFastl + 0.6*rough5ExFastl)
    wT = (1 - np.exp(-0.29*loud5ExZwicker))*(1 - np.exp(-5.49*tonal5ExAures))
    wI = 0.075*impuls5ExHEAD/loud5ExZwicker**-1.334

    sqm = np.maximum(0, 103.08 + 339.49*(wS**2) + 121.88*(wMod**2) + 77.2*(wT**2) + 29.29*(wI**2))

    psych_annoy = loud5ExZwicker*(1 + np.sqrt(sqm))

    return psych_annoy