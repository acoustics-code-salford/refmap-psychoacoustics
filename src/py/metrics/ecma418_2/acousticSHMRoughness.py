# -*- coding: utf-8 -*-
# %% Preamble
"""
acousticSHMRoughness.py
----------------------

Returns roughness values according to ECMA-418-2:2025 using the Sottek Hearing
Model) for an input calibrated single mono or single stereo audio (sound
pressure) time-series signal, p.

Requirements
------------
numpy
scipy
matplotlib
tqdm
bottleneck
acoustic-toolbox
refmap-psychoacoustics (metrics.ecma418_2, dsp.filterFuncs and
                        utils.formatFuncs)

Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford

Date created: 29/05/2023
Date last modified: 01/07/2025
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

"""

# %% Import block
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.fft import (fft, ifft)
from scipy.signal import (hilbert, windows, find_peaks)
from src.py.metrics.ecma418_2.acousticSHMSubs import (shmDimensional,
                                                      shmResample,
                                                      shmPreProc,
                                                      shmOutMidEarFilter,
                                                      shmAuditoryFiltBank,
                                                      shmSignalSegmentBlocks,
                                                      shmSignalSegment,
                                                      shmBasisLoudness,
                                                      shmDownsample,
                                                      shmRoughWeight,
                                                      shmRoughLowPass)
from tqdm import tqdm
import bottleneck as bn
from src.py.dsp.filterFuncs import A_weight_T
from src.py.utils.formatFuncs import roundTrad
from acoustic_toolbox.signal import rms

# %% Module settings
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['mathtext.fontset'] = 'stixsans'

plt.rc('font', size=16)  # controls default text sizes
plt.rc('axes', titlesize=22,
       labelsize=22)  # fontsize of the axes title and x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the tick labels
plt.rc('ytick', labelsize=20)  # fontsize of the tick labels
plt.rc('legend', fontsize=16)  # legend fontsize
plt.rc('figure', titlesize=24)  # fontsize of the figure title


# %% acousticSHMRoughness
def acousticSHMRoughness(p, sampleRateIn, axisN=0, soundField='freeFrontal',
                         waitBar=True, outPlot=False, binaural=True):
    """
    Inputs
    ------
    p : 1D or 2D array
        the input signal as single mono or stereo audio (sound
        pressure) signals

    sampleRateIn : integer
                   the sample rate (frequency) of the input signal(s)

    axisN : integer (0 or 1, default: 0)
            the time axis along which to calculate the roughness

    soundField : keyword string (default: 'freeFrontal')
                 determines whether the 'freeFrontal' or 'diffuse' field stages
                 are applied in the outer-middle ear filter, or 'noOuter' uses
                 only the middle ear stage, or 'noEar' omits ear filtering.
                 Note: these last two options are beyond the scope of the
                 standard, but may be useful if recordings made using
                 artificial outer/middle ear are to be processed using the
                 specific recorded responses.

    waitBar : keyword string (default: True)
              determines whether a progress bar displays during processing
              (set waitBar to false for doing multi-file parallel calculations)

    outPlot : Boolean (default: False)
              flag indicating whether to generate a figure from the output
              (set outplot to false for doing multi-file parallel calculations)

    binaural : Boolean (default: True)
               flag indicating whether to output combined binaural roughness
               for stereo input signal
    Returns
    -------
    roughnessSHM : dict
                  contains the output

    roughnessSHM contains the following outputs:

    specRoughness : 2D or 3D array
                    time-dependent specific roughness for each critical band
                    arranged as [time, bands(, channels)]

    specRoughnessAvg : 1D or 2D array
                      time-averaged specific roughness for each critical band
                      arranged as [bands(, channels)]

    roughnessTDep : 1D or 2D array
                   time-dependent overall roughness
                   arranged as [time(, channels)]

    roughness90pc : 1D or 2D array
                    time-aggregated (90th percentile) overall roughness
                    arranged as [roughness(, channels)]

    bandCentreFreqs : 1D array
                      centre frequencies corresponding with each critical band
                      rate

    timeOut : 1D array
              time (seconds) corresponding with time-dependent outputs

    soundField : string
                 identifies the soundfield type applied (the input argument
                 soundField)

    If outPlot=True, a set of plots is returned illustrating the energy
    time-averaged A-weighted sound level, the time-dependent specific and
    overall roughness, with the latter also indicating the time-aggregated
    value. A set of plots is returned for each input channel.

    If binaural=true, a corresponding set of outputs for the binaural
    loudness are also contained in roughnessSHM.

    Assumptions
    -----------
    The input signal is calibrated to units of acoustic pressure in Pascals
    (Pa).

    Checked by:
    Date last checked:

    """
    # %% Input checks
    # Orient input matrix
    if axisN in [0, 1]:
        if axisN == 1:
            p = p.T
    else:
        raise ValueError("Input axisN must be an integer 0 or 1")

    # Check the length of the input data (must be longer than 300 ms)
    if p.shape[0] <= 300/1000*sampleRateIn:
        raise ValueError('Input signal is too short along the specified axis to calculate tonality (must be longer than 300 ms)')

    # Check the channel number of the input data
    if p.shape[1] > 2:
        raise ValueError('Input signal comprises more than two channels')

    chansIn = p.shape[1]
    if chansIn > 1:
        chans = ["Stereo left",
                 "Stereo right"]
    else:
        chans = ["Mono"]
    # end of if branch for channel number check

    # %% Define constants

    signalT = p.shape[0]/sampleRateIn  # duration of input signal
    # Signal sample rate prescribed to be 48kHz (to be used for resampling), Section 5.1.1 ECMA-418-2:2025 [r_s]
    sampleRate48k = 48e3
    # defined in Section 5.1.4.1 ECMA-418-2:2025 [deltaf(f=0)]
    deltaFreq0 = 81.9289
    # Half-overlapping Bark band centre-frequency denominator constant defined in Section 5.1.4.1 ECMA-418-2:2025
    c = 0.1618

    dz = 0.5  # critical band overlap
    # half-overlapping critical band rate scale [z]
    halfBark = np.arange(0.5, 27, dz)
    nBands = len(halfBark)  # number of critical bands
    # Section 5.1.4.1 Equation 9 ECMA-418-2:2025 [F(z)]
    bandCentreFreqs = (deltaFreq0/c)*np.sinh(c*halfBark)

    # Block and hop sizes Section 6.2.2 Table 4 ECMA-418-2:2025
    overlap = 0.75  # block overlap proportion
    # block size [s_b(z)]
    blockSize = 16384
    # hop sizes (section 5.1.2 footnote 3 ECMA 418-2:2022) [s_h(z)]
    hopSize = int(((1 - overlap)*blockSize))

    # Downsampled block and hop sizes Section 7.1.2 ECMA-418-2:2025
    downSample = 32  # downsampling factor
    sampleRate1500 = sampleRate48k/downSample
    blockSize1500 = int(blockSize/downSample)
    # hopSize1500 = (1 - overlap)*blockSize1500
    resDFT1500 = sampleRate1500/blockSize1500  # DFT resolution (section 7.1.5.1) [deltaf]

    # Modulation rate error correction values Table 8, Section 7.1.5.1
    # ECMA-418-2:2025 [E(theta)]
    errorCorrection = np.array([0.0000, 0.0457, 0.0907, 0.1346, 0.1765, 0.2157, 0.2515,
                                0.2828, 0.3084, 0.3269, 0.3364, 0.3348, 0.3188, 0.2844,
                                0.2259, 0.1351, 0.0000])
    errorCorrection = np.hstack((errorCorrection, np.flip(-errorCorrection[0:-1]), 0))
    
    # High modulation rate roughness perceptual scaling function
    # (section 7.1.5.2 ECMA-418-2:2025)
    # Table 11 ECMA-418-2:2025 [r_1; r_2]
    roughScaleParams = np.vstack(([0.3560, 0.8024], [0.8049, 0.9333]))
    roughScaleParams = np.vstack((roughScaleParams[:, 0]*np.ones([np.sum(bandCentreFreqs < 1e3), 2]),
                                  roughScaleParams[:, 1]*np.ones([np.sum(bandCentreFreqs >= 1e3), 2]))).T
    # Equation 84 ECMA-418-2:2025 [r_max(z)]
    roughScale = 1/(1 + roughScaleParams[0, :]
                    * np.abs(np.log2(bandCentreFreqs/1000))
                    ** roughScaleParams[1, :])
    # Note: this is to ease parallelised calculations
    roughScale = shmDimensional(roughScale, targetDim=3, where='first')

    # High/low modulation rate roughness perceptual weighting function parameters
    # (section 7.1.5.2 ECMA-418-2:2025)
    # Equation 86 ECMA-418-2:2025 [f_max(z)]
    modfreqMaxWeight = 72.6937*(1 - 1.1739*np.exp(-5.4583*bandCentreFreqs/1000))

    # Equation 87 ECMA-418-2:2025 [q_1; q_2(z)]
    roughHiWeightParams = np.vstack((1.2822*np.ones(bandCentreFreqs.size),
                                     0.2471*np.ones(bandCentreFreqs.size)))
    mask = bandCentreFreqs/1000 >= 2**-3.4253
    roughHiWeightParams[1, mask] = 0.2471 + 0.0129*(np.log2(bandCentreFreqs[mask]/1000) + 3.4253)**2
    # Note: this is to ease parallelised calculations
    roughHiWeightParams = np.expand_dims(roughHiWeightParams, axis=1)

    # %% Signal processing

    # Input pre-processing
    # --------------------
    if sampleRateIn != sampleRate48k:  # Resample signal
        p_re, _ = shmResample(p, sampleRateIn)
    else:  # don't resample
        p_re = p
    # end of if branch for resampling

    # Input signal samples
    n_samples = p_re.shape[0]

    # Section 5.1.2 ECMA-418-2:2025 Fade in weighting and zero-padding
    pn = shmPreProc(p_re, blockSize=blockSize, hopSize=hopSize, padStart=True,
                    padEnd=False)

    # Apply outer & middle ear filter
    # -------------------------------
    #
    # Section 5.1.3.2 ECMA-418-2:2025 Outer and middle/inner ear signal filtering
    pn_om = shmOutMidEarFilter(pn, soundField)

    # Loop through channels in file
    # -----------------------------
    if waitBar:
        chanIter = tqdm(range(chansIn), desc="Channels")
    else:
        chanIter = range(chansIn)

    for chan in chanIter:
        # Apply auditory filter bank
        # --------------------------
        # Filter equalised signal using 53 1/2Bark ERB filters according to
        # Section 5.1.4.2 ECMA-418-2:2025
        pn_omz = shmAuditoryFiltBank(pn_om[:, chan])

        # Note: At this stage, typical computer RAM limits impose a need to
        # loop through the critical bands rather than continue with a
        # parallelised approach, until later downsampling is applied

        if waitBar:
            envIter = tqdm(range(nBands), desc="Envelope extraction")
        else:
            envIter = range(nBands)

        # pre-allocate output arrays
        i_start = 0  # index to start segmentation block processing from
        _, nBlocks, _ = shmSignalSegmentBlocks(pn_omz[:, 0],
                                                    blockSize=blockSize,
                                                    overlap=overlap, axisN=0,
                                                    i_start=i_start,
                                                    endShrink=True)
        basisLoudness = np.zeros([nBlocks, nBands])
        envelopes = np.zeros([blockSize1500, nBlocks, nBands])
        for zBand in envIter:
            # Segmentation into blocks
            # ------------------------
    
            # Section 5.1.5 ECMA-418-2:2025
            pn_lz, iBlocksOut = shmSignalSegment(pn_omz[:, zBand],
                                                 blockSize=blockSize,
                                                 overlap=overlap, axisN=0,
                                                 i_start=i_start,
                                                 endShrink=True)

            # Transformation into Loudness
            # ----------------------------
            # Sections 5.1.6 to 5.1.9 ECMA-418-2:2025
            _, bandBasisLoudness, _ = shmBasisLoudness(signalSegmented=pn_lz.copy(),
                                                       bandCentreFreq=bandCentreFreqs[zBand])
            basisLoudness[:, zBand] = bandBasisLoudness

            # Envelope power spectral analysis
            # --------------------------------
            # Sections 7.1.2 ECMA-418-2:2025
            # magnitude of Hilbert transform with downsample - Equation 65
            # [p(ntilde)_E,l,z]
            envelopes[:, :, zBand] = shmDownsample(np.abs(hilbert(pn_lz,
                                                                  axis=0)),
                                                   axisN=0, downSample=32)

        # end of for loop for obtaining low frequency signal envelopes

        # Note: With downsampled envelope signals, parallelised approach can continue

        # Section 7.1.3 equation 66 ECMA-418-2:2025 [Phi(k)_E,l,z]
        modSpectra = np.zeros(envelopes.shape)
        envelopeWin = envelopes*np.tile(shmDimensional(windows.hann(blockSize1500,
                                                                    sym=False),
                                                       targetDim=3),
                                        reps=[1, envelopes.shape[1],
                                              nBands])/np.sqrt(0.375)
        # Equation 66 & 67
        denom = shmDimensional(np.max(basisLoudness, 1))*np.sum(envelopeWin**2,
                                                                0)
        mask = denom != 0  # Equation 66 criteria for masking
        maskRep = np.broadcast_to(mask, modSpectra.shape)  # broadcast mask
        scaling = np.divide(basisLoudness**2, denom, out=np.zeros_like(denom),
                            where=denom != 0)  # Equation 66 factor
        scalingRep = np.broadcast_to(scaling, maskRep.shape)  # broadcast scaling
        modSpectra[maskRep] = np.ravel((np.reshape(scalingRep[maskRep],
                                                   (blockSize1500, -1, nBands))
                                        * np.abs(fft(np.reshape(envelopeWin[maskRep],
                                                                (blockSize1500,
                                                                 -1, nBands)),
                                                     axis=0))**2))

        # Envelope noise reduction
        # ------------------------
        # section 7.1.4 ECMA-418-2:2025
        #TODO change this to numpy slide_tricks
        modSpectraAvg = bn.move_mean(modSpectra, window=3, axis=2)
        modSpectraAvg = np.concatenate((shmDimensional(modSpectra[:, :, 0],
                                                       targetDim=3,
                                                       where='last'),
                                        modSpectraAvg[:, :, 2:],
                                        shmDimensional(modSpectra[:, :, -1],
                                        targetDim=3, where='last')), axis=2)

        modSpectraAvgSum = np.sum(modSpectraAvg, axis=2)  # Equation 68 [s(l,k)]

        # Equation 71 ECMA-418-2:2025 [wtilde(l,k)]
        clipWeight = (0.0856*modSpectraAvgSum[0:int(modSpectraAvg.shape[0]/2) + 1, :]
                      / (np.median(modSpectraAvgSum[2:int(modSpectraAvg.shape[0]/2), :],
                                   axis=0) + 1e-10)
                      * shmDimensional(np.minimum(np.maximum(0.1891*np.exp(0.012*np.arange(0, int(modSpectraAvg.shape[0]/2)
                                                                                           + 1)),
                                                             0),
                                                  1)))

        # Equation 70 ECMA-418-2:2025 [w(l,k)]
        weightingFactor1 = np.zeros(modSpectraAvgSum[0:257, :].shape)
        mask = clipWeight >= 0.05*np.max(clipWeight[2:256, :], axis=0)
        weightingFactor1[mask] = np.minimum(np.maximum(clipWeight[mask] - 0.1407, 0), 1)
        weightingFactor = np.concatenate((weightingFactor1,
                                          np.flipud(weightingFactor1[1:256,
                                                                     :])),
                                         axis=0)

        # Calculate noise-reduced, scaled, weighted modulation power spectra
        # Equation 69 [Phihat(k)_E,l,z]
        modWeightSpectraAvg = modSpectraAvg*shmDimensional(weightingFactor,
                                                           targetDim=3)

        # Spectral weighting
        # ------------------
        # Section 7.1.5 ECMA-418-2:2025
        # theta used in equation 79, including additional index for
        # errorCorrection terms from table 10
        theta = np.arange(34)
        modAmp = np.zeros([10, nBlocks, nBands])
        modRate = np.zeros([10, nBlocks, nBands])

        if waitBar:
            rateIter = tqdm(range(nBands),
                            desc="Critical band modulation rates")
        else:
            rateIter = range(nBands)

        if waitBar:
            blockIter = tqdm(range(nBlocks),
                             desc="Time block modulation rates")
        else:
            blockIter = range(nBlocks)

        for zBand in rateIter:
            # Section 7.1.5.1 ECMA-418-2:2025
            for lBlock in blockIter:
                # identify peaks in each block (for each band)
                startIdx = 2
                endIdx = 255
                modWeightSpectraAvgBandBlock = modWeightSpectraAvg[startIdx:endIdx,
                                                                   lBlock,
                                                                   zBand]
                kLocs, pkProps = find_peaks(modWeightSpectraAvgBandBlock,
                                            prominence=0)
                PhiPks = modWeightSpectraAvgBandBlock[kLocs]
                proms = pkProps['prominences']
    
                # reindex kLocs to match spectral start index used in findpeaks
                # for indexing into modulation spectra matrices
                kLocs = kLocs + startIdx
                
                # we canonly have peaks at k = 3:254
                mask = np.isin(kLocs, range(3, 255))
                kLocs = kLocs[mask]
                PhiPks = PhiPks[mask]
    
                # consider 10 highest prominence peaks only
                if len(proms) > 10:
                    promsSorted = np.sort(proms)[::-1]
                    iiSort = np.argsort(proms)[::-1]
                    mask = proms >= promsSorted[9]
                    
                    # if branch to deal with duplicated peak prominences
                    if sum(mask) > 10:
                       mask = mask[iiSort <= 9]
                    # end of if branch for duplicated peak prominences
    
                    PhiPks = PhiPks[mask]
                    kLocs = kLocs[mask]
    
                # end of if branch to select 10 highest prominence peaks
    
                # consider peaks meeting criterion
                if PhiPks.size != 0:
                    mask = PhiPks > 0.05*np.max(PhiPks)  # Equation 72 criterion
                    PhiPks = PhiPks[mask]  # [Phihat(k_p,i(l,z))]
                    kLocs = kLocs[mask]
                    # loop over peaks to obtain modulation rates
                    
                    if waitBar:
                        peakIter = tqdm(range(len(PhiPks)))
                    else:
                        peakIter = range(len(PhiPks))
                    
                    for iPeak in peakIter:
                        # Equation 74 ECMA-418-2:2025
                        # [Phihat_E,l,z]
                        modAmpMat = np.vstack((modWeightSpectraAvg[kLocs[iPeak] - 1, lBlock, zBand],
                                               modWeightSpectraAvg[kLocs[iPeak], lBlock, zBand],
                                               modWeightSpectraAvg[kLocs[iPeak] + 1, lBlock, zBand]))
                        
                        # Equation 82 [A_i(l,z)]
                        modAmp[iPeak, lBlock, zBand] = np.sum(modAmpMat)
    
                        # Equation 75 ECMA-418-2:2025
                        # [K]
                        modIndexMat = np.vstack((np.hstack(((kLocs[iPeak] - 1)**2, kLocs[iPeak] - 1, 1)),
                                                 np.hstack(((kLocs[iPeak])**2, kLocs[iPeak], 1)),
                                                 np.hstack(((kLocs[iPeak] + 1)**2, kLocs[iPeak] + 1, 1))))
    
                        # Equation 73 solution [C]
                        coeffVec = np.linalg.solve(modIndexMat, modAmpMat)
    
                        # Equation 76 ECMA-418-2:2025 [ftilde_p,i(l,z)]
                        modRateEst = (-(coeffVec[1]/(2*coeffVec[0]))*resDFT1500).item(0)
    
                        # Equation 79 ECMA-418-2:2025 [beta(theta)]
                        errorBeta = ((np.floor(modRateEst/resDFT1500)
                                      + theta[:33]/32)*resDFT1500
                                      - (modRateEst + errorCorrection[theta[:33]]))
    
                        # Equation 80 ECMA-418-2:2025 [theta_min]
                        thetaMinError = np.argmin(np.abs(errorBeta))
    
                        # Equation 81 ECMA-418-2:2025 [theta_corr]
                        if (thetaMinError > 0) and (errorBeta[thetaMinError]*errorBeta[thetaMinError - 1] < 0):
                            thetaCorr = thetaMinError
                        else:
                            thetaCorr = thetaMinError + 1
                        # end of eq 81 if-branch
    
                        # Equation 78 ECMA-418-2:2025
                        # [rho(ftilde_p,i(l,z))]
                        biasAdjust = (errorCorrection[thetaCorr - 1]
                                      - (errorCorrection[thetaCorr]
                                         - errorCorrection[thetaCorr - 1])
                                         * errorBeta[thetaCorr - 1]
                                         / (errorBeta[thetaCorr]
                                            - errorBeta[thetaCorr - 1]))
    
                        # Equation 77 ECMA-418-2:2025 [f_p,i(l,z)]
                        modRate[iPeak, lBlock, zBand] = modRateEst + biasAdjust
    
                    # end of for loop over peaks in block per band
                # end of if branch for detected peaks in modulation spectrum
            # end of for loop over blocks for peak detection      
        # end  of for loop over bands for modulation spectral weighting

        # Section 7.1.5.2 ECMA-418-2:2025 - Weighting for high modulation rates
        # Equation 85 [G_l,z,i(f_p,i(l,z))]
        roughHiWeight = shmRoughWeight(modRate, modfreqMaxWeight,
                                       roughHiWeightParams)

    # Equation 83 [Atilde_i(l,z)]
    modAmpHiWeight = modAmp*roughScale
    mask = modRate <= resDFT1500
    modAmpHiWeight[mask] = 0
    mask = modRate > modfreqMaxWeight
    modAmpHiWeight[mask] = modAmpHiWeight[mask]*roughHiWeight[mask]


    # %% Output plotting

    # Plot figures
    # ------------
    if outPlot:
        if waitBar:
            plotIter = tqdm(range(chansOut), desc="Channels")
        else:
            plotIter = range(chansOut)

        for chan in plotIter:
            # Plot results
            cmap_inferno = mpl.colormaps['inferno']
            chan_lab = chans[chan]
            fig, axs = plt.subplots(nrows=2, ncols=1, figsize=[10.5, 7.5],
                                    layout='constrained')

            ax1 = axs[0]
            pmesh = ax1.pcolormesh(timeOut, bandCentreFreqs,
                                   np.swapaxes(specLoudness[:, :, chan], 0, 1),
                                   cmap=cmap_inferno,
                                   vmin=0,
                                   vmax=np.ceil(np.max(loudnessTDep[:,
                                                                    chan])*10)/10,
                                   shading='gouraud')
            ax1.set(xlim=[timeOut[1],
                          timeOut[-1] + (timeOut[1] - timeOut[0])],
                    xlabel="Time, s",
                    ylim=[bandCentreFreqs[0], bandCentreFreqs[-1]],
                    yscale='log',
                    yticks=[63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],
                    yticklabels=["63", "125", "250", "500", "1k", "2k", "4k",
                                 "8k", "16k"],
                    ylabel="Frequency, Hz")
            ax1.minorticks_off()
            cbax = ax1.inset_axes([1.05, 0, 0.05, 1])
            fig.colorbar(pmesh, ax=ax1,
                         label=(r"Specific loudness,"
                                "\n"
                                r"$\mathregular{tu_{SHM}}/\mathregular{Bark_{SHM}}$"),
                         aspect=10, cax=cbax)

            ax2 = axs[1]
            ax2.plot(timeOut, loudnessPowAvg[chan]*np.ones(timeOut.size),
                     color=cmap_inferno(33/255), linewidth=1,
                     label=("Time-" + "\n" + "average"))
            ax2.plot(timeOut, loudnessTDep[:, chan],
                     color=cmap_inferno(165/255),
                     linewidth=0.75, label=("Time-" + "\n" + "dependent"))
            ax2.set(xlim=[timeOut[0], timeOut[-1] + timeOut[1] - timeOut[0]],
                    xlabel="Time, s",
                    ylim=[0, 1.1*np.ceil(np.max(loudnessTDep[:, chan])*10)/10],
                    ylabel=(r"Loudness, $\mathregular{sone_{SHM}}$"))
            ax2.grid(alpha=0.075, linestyle='--')
            ax2.legend(bbox_to_anchor=(1, 0.85), title="Overall")

            # Filter signal to determine A-weighted time-averaged level
            if chan == 3:
                pA = A_weight_T(p_re, fs=sampleRate48k)
                LAeq2 = 20*np.log10(np.array([rms(pA[:, 0]),
                                              rms(pA[:, 1])])/2e-5)
                # take the higher channel level as representative (PD ISO/TS
                # 12913-3:2019 Annex D)
                LAeq = np.max(LAeq2)
                lr = np.argmax(LAeq2)
                # identify which channel is higher
                if lr == 0:
                    whichEar = " left ear"
                else:
                    whichEar = " right ear"
                # end of if branch to identify which channel is higher

                chan_lab = chan_lab + whichEar
            else:
                pA = A_weight_T(p_re[:, chan], fs=sampleRate48k)
                LAeq = 20*np.log10(rms(pA)/2e-5)

            fig.suptitle(t=(chan_lab + " signal sound pressure level = " +
                            str(roundTrad(LAeq, 1)) +
                            r"dB $\mathregular{\mathit{L}_{Aeq}}$"))
            fig.show()
        # end of for loop over channels
    # end of if branch for plotting

    # %% Output assignment

    # Assign outputs to structure
    if chansOut == 3:
        loudnessSHM = dict()
        loudnessSHM.update({'specLoudness': specLoudness[:, :, 0:2]})
        loudnessSHM.update({'specTonalLoudness': specTonalLoudness[:, :, 0:2]})
        loudnessSHM.update({'specNoiseLoudness': specNoiseLoudness[:, :, 0:2]})
        loudnessSHM.update({'specLoudnessPowAvg': specLoudnessPowAvg[:, 0:2]})
        loudnessSHM.update({'loudnessTDep': loudnessTDep[:, 0:2]})
        loudnessSHM.update({'loudnessPowAvg': loudnessPowAvg[0:2]})
        loudnessSHM.update({'specLoudnessBin': specLoudness[:, :, 3]})
        loudnessSHM.update({'specTonalLoudnessBin': specTonalLoudness[:, :, 3]})
        loudnessSHM.update({'specNoiseLoudnessBin': specNoiseLoudness[:, :, 3]})
        loudnessSHM.update({'specLoudnessPowAvgBin': specLoudnessPowAvg[:, 3]})
        loudnessSHM.update({'loudnessTDepBin': loudnessTDep[:, 3]})
        loudnessSHM.update({'loudnessPowAvgBin': loudnessPowAvg[3]})
        loudnessSHM.update({'bandCentreFreqs': bandCentreFreqs})
        loudnessSHM.update({'timeOut': timeOut})
        loudnessSHM.update({'soundField': soundField})
    else:
        loudnessSHM.update({'specLoudness': specLoudness})
        loudnessSHM.update({'specTonalLoudness': specTonalLoudness})
        loudnessSHM.update({'specNoiseLoudness': specNoiseLoudness})
        loudnessSHM.update({'specLoudnessPowAvg': specLoudnessPowAvg})
        loudnessSHM.update({'loudnessTDep': loudnessTDep})
        loudnessSHM.update({'loudnessPowAvg': loudnessPowAvg})
        loudnessSHM.update({'bandCentreFreqs': bandCentreFreqs})
        loudnessSHM.update({'timeOut': timeOut})
        loudnessSHM.update({'soundField': soundField})

# end of acousticSHRoughness function