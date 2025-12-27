"""
detectability.py
----------------

Psychoacoustic detectability analysis for audio signals:

psych_detect:
    Returns detectability and discounted sound levels from input target
    source and masker sound pressure time-series signals based on the detectability
    model originally developed at Bolt, Beranek and Newman consulting engineers,
    and developed further at NASA, with discounting of target source sound levels
    using the technique developed by NASA (see References).

Requirements
------------
numpy
scipy
bottleneck

References
----------    
Fidell, S et al, 1974. Prediction of aural detectability of noise
signals, Human Factors 16(4), 373-383.
https://doi.org/10.1177/001872087401600405

Sneddon, M et al, 2003. Laboratory study of the noticeability and
annoyance of low signal-to-noise ratio sounds, Noise Control Engineering
Journal, 51(5), 300-305.
https://doi.org/10.3397/1.2839726

Christian, A, 2021. A construction for the prediction of noise-induced
annoyance in the presence of auditory masking, 181st Meeting of the
Acoustical Society of America.
https://ntrs.nasa.gov/citations/20210024824

Rizzi, SA et al, 2024. Annoyance model assessments of urban air mobility
vehicle operations. 30th AIAA/CEAS Aeroacoustics Conference, Rome, Italy,
June 4-7, 2024. https://doi.org/10.2514/6.2024-3014


Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford

Date created: 05/11/2024
Date last modified: 16/12/2025
Python version: 3.11

Copyright statement: This file and code is part of work undertaken within
the RefMap project (www.refmap.eu), and is subject to licence as detailed
in the code repository
(https://github.com/acoustics-code-salford/refmap-psychoacoustics)

As per the licensing information, please be aware that this code is WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

"""
import numpy as np
import warnings
from math import gcd
import bottleneck as bn
from scipy.signal import butter, resample_poly, sosfiltfilt, sosfilt


# psych_detect function
def psych_detect(signal_target, samp_rate_target, axis_target,
                 signal_masker, samp_rate_masker, axis_masker,
                 time_skip=[0, 0], time_step=0.5, noct=3,
                 freq_range=[20, 20000]):
    """
    Returns detectability and discounted sound levels from input target
    source and masker sound pressure time-series signals.
    
    Parameters
    ----------
    signal_target : ndarray
        Target signal as single mono or stereo audio (sound pressure) signals
        [time, channels].

    sample_rate_target : int
        Sampling frequency of target signal (Hz).
    
    axis_target : int
        Time axis of target signal (0 = rows, 1 = columns).
    
    signal_masker : ndarray
        Masker signal [time, channels].
    
    samp_rate_masker : int
        Sampling frequency of masker signal (Hz).
    
    axis_masker : int
        Time axis of masker signal (0 = rows, 1 = columns).

    time_skip : vector (default: [0, 0])
        Time (seconds) to skip from input signals for calculating
        time-aggregated outputs. [start_skip, end_skip] ignores
        start_skip seconds of the start, and end_skip seconds of the end.

    time_step : number (default: 0.5)
        Time window (seconds) to use for calculating target
        detectability.

    noct : integer (default: 3)
        Number of fractional-octave bands to use (e.g., default 3 = 1/3-octave).
        Note that fractional-octave bands other than 1/3-octave have not been
        validated.

    freq_range : vector (default: [20, 20000])
        The frequency range over which to determine
        detection and discounted spectra (1/3-octave band
        centre-frequencies within this range will be included). [f_min, f_max]

    Returns
    -------
    detectability : dict
        Contains the outputs.
    
    detectability contains the following outputs:

    lpz_spec_target : matrix
        Sound pressure level spectrogram for the input target
        signal, with dimensions [time_out, freq_bands, target_chans].

    lpz_spec_masker : matrix
        Sound pressure level spectrogram for the input masker
        signal, with dimensions [time_out, freq_bands, masker_chans]. 

    lpz_spec_disc_target : matrix
        Sound pressure level spectrogram for the input target
        signal, detectability-discounted, with dimensions
        [time_out, freq_bands, target_chans].

    lpa_tdep_target : matrix or vector
        Time-dependent A-weighted sound pressure level for the
        input target signal, with dimensions [time_out, target_chans].

    lpa_tdep_masker : matrix or vector
        Time-dependent A-weighted sound pressure level for the
        input masker signal, with dimensions [time_out, target_chans].

    lpa_tdep_disc_target : matrix or vector
        Time-dependent A-weighted sound pressure level for the
        input target signal detectability-discounted, with
        dimensions [time_out, target_chans].

    lpa_tdep_discount : vector
        Time-dependent detectability discount values for the
        A-weighted sound pressure levels of target signal vs
        masker signal, with dimensions [time_out, target_chans].

    lae_target : vector
        Overall A-weighted sound exposure level for each input target
        signal channel.

    lae_masker : vector
        Overall A-weighted sound exposure level for each input masker
        signal channel.

    lae_disc_target : vector
        Overall detectability-discounted A-weighted sound
        exposure level for each input target signal channel.

    laeq_target : vector
        Overall A-weighted time-averaged sound level for each input
        target signal channel.

    laeq_masker : vector
        Overall A-weighted time-averaged sound level for each input
        masker signal channel.

    laeq_disc_target : vector
        Overall detectability-discounted A-weighted
        time-averaged sound level for each input target signal
        channel.

    dba_discount : vector
        Overall detectability discount for A-weighted levels for
        each input target signal channel.

    detectability : matrix
        Detectability spectrogram for the input target signal
        vs masker signal, with dimensions [time_out, freq_bands, target_chans]

    detect_tdep_max : matrix or vector
        band-maximum time-dependent detectability, with
        dimensions [time_out, target_chans]

    detect_tdep_int : matrix or vector
        band-integrated time-dependent detectability, with
        dimensions [time_out, target_chans]

    detect_max : vector
        overall maximum detectability for each input target
        signal channel

    detect_int : vector
        overall integrated detectability for each input target
        signal channel

    freq_bands : vector
        1/3-octave band centre-frequencies for input freq_range

    time_out : vector
        the window centre times for the spectrograms.
 
    Assumptions
    -----------
    Both input signals are calibrated to absolute acoustic pressure (Pascals),
    or, otherwise, equal scaling.

    """ 
    if axis_target == 1:
        signal_target = np.transpose(signal_target)
    
    if axis_masker == 1:
        signal_masker = np.transpose(signal_masker)

    # check input signal sampling frequencies are equal, otherwise resample to
    # match higher rate
    if samp_rate_masker > samp_rate_target:
        # upsampled sampling frequency
        samp_rate = samp_rate_masker
        up = int(samp_rate/gcd(samp_rate, samp_rate_target))
        down = int(samp_rate_target/gcd(samp_rate, samp_rate_target))
        signal_target = signal.resample_poly(signal_target, up, down, padtype='line')
        
    elif samp_rate_target > samp_rate_masker:
        samp_rate = samp_rate_target
        up = int(samp_rate/gcd(samp_rate, samp_rate_masker))
        down = int(samp_rate_masker/gcd(samp_rate, samp_rate_masker))
        signal_masker = signal.resample_poly(signal_masker, up, down, padtype='line')
    else:
        samp_rate = samp_rate_target
  
    time_steps = int(time_step*samp_rate)  # time step in samples
    
    # check signal lengths match
    if signal_target.shape[0] != signal_masker.shape[0]:
        raise ValueError("The lengths of the input signals on the time axes do not match.")
    
    # check number of channels match, if not, raise a warning and try to broadcast masker
    # to match target channels
    if signal_target.shape[1] != signal_masker.shape[1]:
        warnings.warn("The numbers of channels in the input signals do not match. Attempting to broadcast masker to match target channels (if this behaviour is not expected, review the input signals).", UserWarning)
        try:
            signal_masker = np.tile(signal_masker, (1, signal_target.shape[1]))
            numChans = signal_target.shape[1]
        except:
            raise ValueError("The numbers of channels in the input signals do not match, and broadcasting failed.")
    else:
        numChans = signal_target.shape[1]

    # if target or masker signal array dimensions are 1D, expand to 2D
    if signal_target.ndim == 1:
        signal_target = np.expand_dims(signal_target, axis=1)
    if signal_masker.ndim == 1:
        signal_masker = np.expand_dims(signal_masker, axis=1)

    # ensure F-contiguous arrays for processing speed
    signal_target = np.asfortranarray(signal_target)
    signal_masker = np.asfortranarray(signal_masker)

    # clamp freq_range to valid range
    if freq_range[0] < 20:
        freq_range[0] = 20
        warnings.warn("\nLower frequency range limited to 20 Hz.", UserWarning)
    if freq_range[1] > samp_rate/2:
        freq_range[1] = samp_rate/2
        warnings.warn("\nUpper frequency range limited to Nyquist frequency.", UserWarning)

    # calculate fractional-octave band LZeq(t)  
    fm, _, _ = noctf(fl=freq_range[0], fh=freq_range[1], n=noct)

    # pre-allocate filtered signal array
    signal_tob = np.zeros((signal_target.shape[0], signal_target.shape[1], fm.size), order='F')
    # filter target and masker signals into fractional-octave band LZeq(t)
    for signal, targ_mask_idx in enumerate([signal_target, signal_masker]):
        for f, f_idx in enumerate(fm):
            signal_tob[:, :, f_idx] = noct_filter(x=signal, n=noct, fm=f,
                                                  fs=samp_rate, order=6,
                                                  axis=0, fwd_bwd=True)
            
            signal_tob_rms = np.sqrt(bn.move_mean(np.square(signal_tob),
                                     window=time_steps, axis=0))

        if targ_mask_idx == 0:
            signal_target_tob_leq = 20*np.log10(signal_tob_rms[time_steps - 1::time_steps, :, :]/2e-5)
        else:
            signal_masker_tob_leq = 20*np.log10(signal_tob_rms[time_steps - 1::time_steps, :, :]/2e-5)


    # calculate detectability and discounted levels
    detectability = psych_detect_from_lzeqt(lzeqt_spec_target=signal_target_tob_leq,
                                            lzeqt_spec_masker=signal_masker_tob_leq,
                                            time_step=time_step,
                                            time_skip=time_skip,
                                            freq_band_range=freq_range)

    return detectability

    
def psych_detect_from_lzeqt(lzeqt_spec_target, lzeqt_spec_masker,
                            axis_target=0, axis_masker=0,
                            time_step=0.5, time_skip=[0, 0],
                            noct=3, freq_band_range=[20, 20000]):
    """
    Returns detectability and discounted sound levels from input target
    source and masker unweighted (Z-weighted) short-term equivalent continuous
    (energy-time-averaged) sound pressure level LZeq(t) 1/3-octave-band
    spectrograms.

    Parameters
    ----------
    lzeqt_spec_target : 2D or 3D array
        Target signal LZeq(t) spectrogram [time, freq_bands, channels].
    
    lzeqt_spec_masker : 2D or 3D array
        Masker signal LZeq(t) spectrogram [time, freq_bands, channels].

    axis_target : int (default: 0)
        Time axis of target signal spectrogram (0 = rows, 1 = columns).
    
    axis_masker : int (default: 0)
        Time axis of masker signal spectrogram (0 = rows, 1 = columns).
    
    time_step : number (default: 0.5)
        Time step (seconds) corresponding with the input LZeq(t) spectrograms.
        This is also assumed to be the averaging time for the LZeq(t) values.
    
    time_skip : 1D array-like (default: [0, 0])
        Time (seconds) to ignore from input signals for calculating
        time-aggregated outputs [start_skip, end_skip].

    noct : integer (default: 3)
        Number of fractional-octave bands to use (e.g., default 3 = 1/3-octave).
        Note that fractional-octave bands other than 1/3-octave have not been
        validated.
    
    freq_band_range : 1D array-like (default: [20, 20000])
        The frequency band range corresponding with the input spectra (must match
        both target and masker) [f_min, f_max].

    Returns
    -------
    detectability : dict
        Contains the outputs.
    
    detectability contains the following outputs:

    lpz_spec_target : matrix
        Sound pressure level spectrogram for the input target
        signal, with dimensions [time_out, freq_bands, target_chans].

    lpz_spec_masker : matrix
        Sound pressure level spectrogram for the input masker
        signal, with dimensions [time_out, freq_bands, masker_chans]. 

    lpz_spec_disc_target : matrix
        Sound pressure level spectrogram for the input target
        signal, detectability-discounted, with dimensions
        [time_out, freq_bands, target_chans].

    lpa_tdep_target : matrix or vector
        Time-dependent A-weighted sound pressure level for the
        input target signal, with dimensions [time_out, target_chans].

    lpa_tdep_masker : matrix or vector
        Time-dependent A-weighted sound pressure level for the
        input masker signal, with dimensions [time_out, target_chans].

    lpa_tdep_disc_target : matrix or vector
        Time-dependent A-weighted sound pressure level for the
        input target signal detectability-discounted, with
        dimensions [time_out, target_chans].

    lpa_tdep_discount : vector
        Time-dependent detectability discount values for the
        A-weighted sound pressure levels of target signal vs
        masker signal, with dimensions [time_out, target_chans].

    lae_target : vector
        Overall A-weighted sound exposure level for each input target
        signal channel.

    lae_masker : vector
        Overall A-weighted sound exposure level for each input masker
        signal channel.

    lae_disc_target : vector
        Overall detectability-discounted A-weighted sound
        exposure level for each input target signal channel.

    laeq_target : vector
        Overall A-weighted time-averaged sound level for each input
        target signal channel.

    laeq_masker : vector
        Overall A-weighted time-averaged sound level for each input
        masker signal channel.

    laeq_disc_target : vector
        Overall detectability-discounted A-weighted
        time-averaged sound level for each input target signal
        channel.

    dba_discount : vector
        Overall detectability discount for A-weighted levels for
        each input target signal channel.

    detectability : matrix
        Detectability spectrogram for the input target signal
        vs masker signal, with dimensions [time_out, freq_bands, target_chans]

    detect_tdep_max : matrix or vector
        band-maximum time-dependent detectability, with
        dimensions [time_out, target_chans]

    detect_tdep_int : matrix or vector
        band-integrated time-dependent detectability, with
        dimensions [time_out, target_chans]

    detect_max : vector
        overall maximum detectability for each input target
        signal channel

    detect_int : vector
        overall integrated detectability for each input target
        signal channel

    freq_bands : vector
        1/3-octave band centre-frequencies for input freq_range

    time_out : vector
        the window centre times for the spectrograms.
 
    Assumptions
    -----------
    Input spectra are assumed to be 1/3-octave bands unweighted LZeq(t).
    The default time_step=0.5 seconds corresponds with the NASA approach.
    Alternative time or frequency resolutions have not been validated.
    """

    # %% Define constants
    efficiency_factor = 0.3  # \eta
    target_detect = 2  # target d'
    discount_half_power = 3  # \alpha, d_b
    detect_knee = 14  # \delta, d' value at which 3 d_b knee occurs
    discount_rate = 1  # \rho, rate at which discount function dimishes with reducing detectability

    # diffuse field narrowband noise hearing thresholds for 18-25 year-olds with
    # normal hearing from ISO 389-7:2019 (20 Hz - 16 kHz)
    # NOTE: watch out for non-standard frequencies in the standard Table 1!
    hear_thresholds_df3897 = np.array([78.1, 68.7, 59.5, 51.1, 44.0, 37.5, 31.5, 26.5, 22.1,
                                       17.9, 14.4, 11.4, 8.4, 5.8, 3.8, 2.1, 1.0, 0.8, 1.9,
                                       0.5, -1.5, -3.1, -4.0, -3.8, -1.8, 2.5, 6.8, 8.4, 14.4,
                                       43.7])

    # estimated full range diffuse field narrowband noise hearing thresholds
    # for 18-25 year-olds with normal hearing 
    hear_thresholds_df = np.append(hear_thresholds_df3897, 70.4)

    # 1/3-octave A-weighting dB values (20 Hz - 20 kHz)
    a_weight = np.array([-50.5, -44.7, -39.4, -34.6, -30.2, -26.2, -22.5, -19.1, -16.1,
                         -13.4, -10.9, -8.6, -6.6, -4.8, -3.2, -1.9, -0.8, 0.0, 0.6, 1.0,
                         1.2, 1.3, 1.2, 1.0, 0.5, -0.1, -1.1, -2.5, -4.3, -6.6, -9.3])

    f, _, _ = noctf(fl=freq_band_range[0], fh=freq_band_range[1], n=noct)

    # convert timeSkip to timeStep indices
    itime_skip = time_skip/time_step

    # Equivalent rectangular bandwidth for auditory filter (Glasberg & Moore,
    # 1990, equation 3)
    erb_bandwidth = 24.7*(1 + 4.37/1000*f)

    # Calculate 1/3 octave bandwidths according to BS EN IEC 61640-1:2014
    g10 = 10**(3/10)  # octave ratio coefficient (base-ten)
    oct_ratio = g10**(0.5/noct)  # octave ratio

    f1 = f/oct_ratio  # output range of exact lower band-edge frequencies
    f2 = f*oct_ratio  # output range of exact upper band-edge frequencies
    f_bandwidth = f2 - f1

    # Calculate high frequency weighting "kick" factor
    w_kick = np.zeros(f.shape)
    w_kick(f >= 2500) = -6*(np.log10(f(f >= 2500)/2500))**2

    # Calculate detection efficiency
    detect_efficiency = efficiency_factor*np.sqrt(time_step)*np.sqrt(f_bandwidth)*np.sqrt(f_bandwidth/erb_bandwidth)*10**(w_kick/10)

    # Calculate equivalent auditory system noise
    ind = range(-17, 13)  # range of frequency indices for hearing threshold bands
    ind_vis = range(-17, 14)  # used for visualisations
    if np.mod(noct, 1) == 0:
        fm = g10**(ind/noct)*1000
        fm_vis = g10**(ind_vis/noct)*1000
    else:
        fm = g10**((2*ind + 1)/(2*noct))*1000
        fm_vis = g10**((2*ind_vis + 1)/(2*noct))*1000

    f_bandwidth_vis = fm_vis*oct_ratio - fm_vis/oct_ratio

    il = np.argmin(np.abs(fm - f[0]))  # find nearest lower exact frequency
    ih = np.argmin(np.abs(fm - f[-1]))  # find nearest higher exact frequency

    eq_auditory_noise = detect_efficiency*(2e-5*10**(hear_thresholds_df[il:ih]/20))**2/target_detect

    # Calculate detectability and discount
    detectability = detect_efficiency*lzeqt_spec_target/(np.tile(lzeqt_spec_masker, (1, 1, rep_masker)) + eq_auditory_noise)
    detect_discount = discount_half_power/(detectability/detect_knee)**discount_rate
    lzeqt_dscnt_spec_target = lzeqt_spec_target - detect_discount

    # Calculate aggregated detectability
    detect_t_max = np.squeeze(np.max(detectability, axis=0))
    detect_t_int = np.squeeze(np.sqrt(np.sum(detectability**2, axis=0)))

    # A-weight time-dependent spectra
    laeq_spec_target = lzeqt_spec_target + a_weight[il:ih]
    laeq_spec_masker = lzeqt_spec_masker + a_weight[il:ih]

    # A-weighted time-dependent spectral power
    pow_a_target = (2e-5*10**(laeq_spec_target/20))**2
    pow_a_masker = (2e-5*10**(laeq_spec_masker/20))**2

    detectability = {}
    detectability['lzeqt_spec_target'] = lzeqt_spec_target
    detectability['lzeqt_spec_masker'] = lzeqt_spec_masker

    return detectability


def noctf(fl, fh, n):
    """
    Return consecutive range of (exact) mid-frequencies, and lower and upper
    band-edge frequencies for 1/n fractional-octave-band filters according to
    BS EN IEC 61260-1:2014

    Parameters
    ----------
    fl : number
      Lowest frequency band of interest.

    fh : number
      Highest frequency band of interest.

    n : integer
      Denominator of the octave fraction 1/n.

    Returns
    -------
    fm : 1D array
      Fractional octave band mid frequencies.

    f1 : 1D array
      Fractional octave band lower band-edge frequencies.

    f2 : 1D array
      Fractional octave band upper band-edge frequencies.

    Assumptions
    -----------

    """
    if not isinstance(n, int):
        raise TypeError("\nOctave fraction denominator must be an integer")

    ind = np.arange(-15*n, 15*n, 1)  # range of frequency indices
    g10 = 10**(3/10)  # octave ratio coefficient (base-ten)
    oct_ratio = g10**(0.5/n)  # octave ratio

    if np.mod(n, 1) == 0:
        f = g10**(ind/n)*1000
    else:
        f = g10**((2*ind + 1)/(2*n))*1000

    il = np.abs(f - fl).argmin()  # find nearest lower exact frequency
    ih = np.abs(f - fh).argmin()  # find nearest higher exact frequency

    fm = f[il:ih+1]  # output range of exact fractional frequencies
    f1 = fm/oct_ratio  # output range of exact lower band-edge frequencies
    f2 = fm*oct_ratio  # output range of exact upper band-edge frequencies

    return fm, f1, f2


def noctf_fbw(f, fbw_order, fs):
    """
    Return noct bandpass filter cutoff frequencies adjusted (pre-warped) to
    compensate for zero-phase forwards-backwards Butterworth filter response

    Developed from filter response functions for (one-way) orders 1-64
    (ie two-way orders 2-128) using linear regression over values adjusted to
    give desired output.

    Parameters
    ----------
    
    f : float or array of floats
        Input frequency which is to be adjusted/warped.

    fbw_order : integer
        Design target order of the fowards-backwards filter, ie 2x the
        value of the order input to the Butterworth filter function.

    fs : integer
        Filter sampling frequency.

    Returns
    -------
    f : float or array of floats
        Adjusted frequencies.

    """
    if not isinstance(fbw_order, int):
        raise TypeError("\nForwards-backwards filter order parameter "
              "'fbw_order' must be an integer")

    if fbw_order < 2 or fbw_order > 128:
        raise ValueError("\nForwards-backwards filter order parameter "
              "'fbw_order' must take a value in the range 2-128")

    K = np.log(fbw_order)

    c = np.array([-3.5676481792823e-05, 1.02690386340958e-03,
                  -1.05748495945522e-02, 5.20134108353335e-02,
                  -1.27978679123245e-01, 1.13201])
    f1_fbw = f[0]/(np.sum(c*K**np.arange(5, -1, -1)))
    f2_fbw = f[1]*(np.sum(c*K**np.arange(5, -1, -1)))
    
    return f1_fbw, f2_fbw


def noctfiltc(n, fm, fs, order=6, fwd_bwd=True):
    """
    Return coefficients in second-order-sections for a 1/n fractional-octave-band
    filter according to BS EN IEC 61260-1:2014.
        
    Parameters
    ----------
    n : integer
        Number defining the octave fraction 1/n.

    fm : float
        Mid-frequency of 1/n octave filter.

    fs : integer
        Filter sampling frequency.

    order : integer
        Filter order.

    fwd_bwd : Boolean (default: True)
        Determines whether the filter is designed as one-way or two-way
        (forwards-backwards, zero-phase). For zero-phase filtering, the
        input parameter order is interpreted as the intended overall filter
        order, so an odd input order is adjusted to be even, with a warning
        (since half the order parameter value is used to define the
        single-pass filter).
 
    Returns
    -------
    sos : array
          the filter coefficients in second-order section form
    zi : float
         the initial condition for the filter

    """
    N = [1, 2, 3, 6, 12, 24, 48]
    if n not in N:
        raise ValueError(("\nOctave fraction denominator must take a value from: "
                          + str(N)))

    # Check and force order to be integer with warning
    if not isinstance(order, int):
        ordernonint = order
        order = int(order)
        warnings.warn(("\nCaution: Butterworth filter must take integer order. Input order "
               + str(ordernonint) + " changed to integer order " + str(order)), UserWarning)

    if fwd_bwd:
        if np.remainder(order, 2) > 0:
            inorder = order
            order = order + 1
        warnings.warn(("\nCaution: for forward-backward processing, even input " \
                   "(target) order is required. Odd input order "
                   + str(inorder) + " changed to even order " + str(order)))

    g10 = 10**(3/10)
    oct_ratio = g10**(0.5/n)

    f1 = fm/oct_ratio
    f2 = fm*oct_ratio

    if fwd_bwd is True:
        # Pre-warping of cutoff frequencies to compensate for forward-backward
        fbw_order = order
        order = int(order/2)
        f1, f2 = noctf_fbw([f1, f2], fbw_order, fs, 'bandpass')

    sos = butter(order, [f1, f2], btype='bandpass', output='sos', fs=fs)
    # initial condition for steady-state of step response to pass to filter function
    zi = sosfilt_zi(sos)

    return sos, zi


def noct_filter(x, n, fm, fs, order=6, axis=0, fwd_bwd=False):
    """
    Return 1/n fractional-octave-band-filtered signals for a single passband
    using second-order-sections.

    Parameters
    ----------
    x : 1D or 2D array
        Input time signals to be filtered.

    n : integer
        Number defining the octave fraction 1/n.

    fm : number
        Mid-frequency for the fractional-octave-band filter.

    fs : number
        Sampling frequency of the input signal.

    order : integer
        Intended order for Butterworth IIR filter design

    axis : integer
        Input signal axis along which to apply filter.

    fwd_bwd : Boolean
        Indicates forward-backward (zero-phase) filtering to be used.
        In this case, 'order' is interpreted as the overall desired target
        filter order, and a warning is raised if 'order' is odd (as the target
        order must be even if forward-backward filtering is used). The filter
        is padded by a constant length taken from the smaller of: i) the length
        of x along 'axis' minus 1, or ii) a value calculated according to [2].

    check : Boolean
        Flag to check for filter compliance with BS EN IEC 61260:2014 [1].

    Returns
    -------
    y : 1D or 2D array
      Filtered output signals with shape: x.shape.


    References
    ----------
    1. BSI, 2014. BS EN IEC 61260-1:2014
    2. Boore, DM, 2005. On pads and filters: processing strong-motion data.
       Bulletin of the Seismological Society of America, 95(2), 745-750.
    """
    N = [1, 2, 3, 6, 12, 24, 48]
    if n not in N:
        raise ValueError(("\nOctave fraction denominator must take a value from: "
                          + str(N)))

    fm, _, f2 = noctf(fm, fm, n)  # convert to exact filter mid-frequency

    # derive filter coefficients
    sos, _ = noctfiltc(n, fm, fs, order=order, fwd_bwd=fwd_bwd)

    if fwd_bwd:
        # filter signal using forward-backward processing
        # generate fbw pad length based on Boore, 2005
        padN = np.min([x.shape[axis]-1, int(x.shape[axis]*1.5*order/2/f2.item())])
        y_filt = sosfiltfilt(sos, x, axis=axis, padtype='constant',
                             padlen=padN)
    else:
        y_filt = sosfilt(sos, x, axis=axis)  # filter signal

    y = y_filt

    return y