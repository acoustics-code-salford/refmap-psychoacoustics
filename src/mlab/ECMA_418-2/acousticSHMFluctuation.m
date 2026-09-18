function fluctuationSHM = acousticSHMFluctuation(p, sampleRateIn, axisN, soundField, waitBar, outPlot, binaural, diagOn)
% fluctuationSHM = acousticSHMFluctuation(p, sampleRateIn, axisN, soundField, waitBar, outPlot, binaural)
%
% Returns fluctuation strength values according to ECMA-418-2:2025
% (using the Sottek Hearing Model) for an input calibrated single mono
% or single stereo audio (sound pressure) time-series signal, p. For stereo
% signals, the binaural fluctuation strength can be calculated, and each
% channel is also analysed separately.
%
% Inputs
% ------
% p : vector or 2D matrix
%   Input signal as single mono or stereo audio (sound
%   pressure) signals
%
% sampleRateIn : integer
%   Sample rate (frequency) of the input signal(s)
%
% axisN : integer (1 or 2, default: 1)
%   Time axis along which to calculate the fluctuation strength
%
% soundField : keyword string (default: 'freeFrontal')
%   Determines whether the 'freeFrontal' or 'diffuse' field stages
%   are applied in the outer-middle ear filter, or 'noOuter' uses
%   only the middle ear stage, or 'noEar' omits ear filtering.
%   Note: these last two options are beyond the scope of the
%   standard, but may be useful if recordings made using
%   artificial outer/middle ear are to be processed using the
%   specific recorded responses.
%
% waitBar : keyword string (default: true)
%   Determines whether a progress bar displays during processing
%
% outPlot : Boolean true/false (default: false)
%   Flag indicating whether to generate a figure from the output
%
% binaural : Boolean true/false (default: true)
%   Flag indicating whether to output combined binaural fluctuation
%   strength for stereo input signal.
%
% diagOn : Boolean true/false (default: false)
%   Diagnostic instrumentation flag, NOT part of the standard's
%   algorithm. When true, one row is logged for every non-quiet
%   block/band processed by the HSA pipeline (Sections 9.1.4-9.1.10),
%   recording the size and conditioning of the linear system solved in
%   Section 9.1.4's initial (potentially many-candidate) fit, together
%   with the block's key intermediate and final values. This is
%   intended purely to help isolate numerical or logical issues (e.g.
%   rank-deficiency of the HSA linear system when many candidate lines
%   are found simultaneously) without needing to re-instrument the code
%   by hand. See fluctuationSHM.diagLog below. Leaving this false (the
%   default) avoids the associated memory/time overhead.
%
% Returns
% -------
%
% fluctuationSHM : structure
%   contains the output
%
% fluctuationSHM contains the following outputs:
%
% specFluctuation : matrix
%   time-dependent specific fluctuation strength for each critical band
%   arranged as [time, bands(, channels)]
%
% specFluctuationAvg : matrix
%   time-averaged specific fluctuation strength for each critical band
%   arranged as [bands(, channels)]
%
% fluctuationTDep : vector or matrix
%   time-dependent overall fluctuation strength arranged as [time(, channels)]
% 
% fluctuation90Pc : number or vector
%   time-aggregated (90th percentile) overall fluctuation strength
%   arranged as [fluctuation strength(, channels)]
%
% bandCentreFreqs : vector
%   centre frequencies corresponding with each critical band rate
%
% timeOut : vector
%   time (seconds) corresponding with time-dependent outputs
%
% soundField : string
%   identifies the soundfield type applied (the input argument
%   soundField)
%
% diagLog : structure (only present if diagOn = true)
%   diagnostic log, not part of the standard's algorithm - see the
%   diagOn input description above. Contains equal-length vectors, one
%   element per logged block/band/channel:
%     chan, zBand, lBlock : indices identifying the block
%     status  : categorical string - 'ok', 'no_modulation' (Section
%               9.1.5, no local maximum and no local minimum),
%               'no_survivors' (Equation 146 threshold left nothing),
%               'no_harmonic_group' (Section 9.1.8 found no valid
%               harmonic grouping), or 'discarded_low_freq' (Section
%               9.1.7, f_c,1,opt < 0.125 Hz)
%     Mc, KL, nUnknown, rcondA, usedPinv : from the initial (potentially
%               many-candidate) HSA fit's diagInfo output (see shmHSA.m)
%     fcOpt   : fine-tuned dominant modulation rate [Hz] (NaN if not
%               reached)
%     nHarm   : number of components in the retained harmonic complex
%               (NaN if not reached)
%     Ahat, powerSum, N_HSA : Section 9.1.9/9.1.10 intermediate values
%               (NaN if not reached)
%     A_lz    : the pre-threshold value of A(l,z) (Equation 159), i.e.
%               before the Section 9.1.10 5.2519 threshold is applied
%               (NaN if not reached)
%     nzb, nze : the Section 9.1.3 envelope analysis window parameters
%               (number of zeros at the start/end of the block) actually
%               used for this block/band's HSA fit
%
% If binaural=true, a corresponding set of outputs for the binaural
% fluctuation strength is also contained in fluctuationSHM
%
% If outplot=true, a set of plots is returned illustrating the energy
% time-averaged A-weighted sound level, the time-dependent specific and
% overall fluctuation strength, with the latter also indicating the
% time-aggregated value. A set of plots is returned for each input channel,
% with another set for the binaural fluctuation strength, if binaural=true.
% In that case, the indicated sound level corresponds with the channel with
% the highest sound level.
%
% Assumptions
% -----------
% The input signal is calibrated to units of acoustic pressure in Pascals
% (Pa).
%
% Requirements
% ------------
% Signal Processing Toolbox
% Audio Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 16/05/2025
% Date last modified: 18/09/2026
% MATLAB version: 2023b
%
% Copyright statement: This file and code is part of work undertaken within
% the RefMap project (www.refmap.eu), and is subject to licence as detailed
% in the code repository
% (https://github.com/acoustics-code-salford/refmap-psychoacoustics)
%
% As per the licensing information, please be aware that this code is
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
% This code calls sub-component file 'cmap_inferno.txt'. The contents of
% the file includes a copy of data obtained from the repository 
% https://github.com/BIDS/colormap, and is CC0 1.0 licensed for modified
% use, see https://creativecommons.org/publicdomain/zero/1.0 for
% information.
%
% Checked by:
% Date last checked:
%
%% Arguments validation
    arguments (Input)
        p (:, :) double {mustBeReal}
        sampleRateIn (1, 1) double {mustBePositive, mustBeInteger}
        axisN (1, 1) {mustBeInteger, mustBeInRange(axisN, 1, 2)} = 1
        soundField (1, :) string {mustBeMember(soundField,...
                                                       {'freeFrontal',...
                                                        'diffuse', ...
                                                        'noOuter', ...
                                                        'noEar'})} = 'freeFrontal'
        waitBar {mustBeNumericOrLogical} = true
        outPlot {mustBeNumericOrLogical} = false
        binaural {mustBeNumericOrLogical} = true
        diagOn {mustBeNumericOrLogical} = false
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Input checks
% Orient input matrix
if axisN == 2
    p = p.';
end

% Check the length of the input data (must be longer than the block used
% for fluctuation strength segmentation, i.e. > 1.3653 s)
if size(p, 1) <= 65536/48e3*sampleRateIn
    error("Error: Input signal is too short along the specified axis to calculate fluctuation strength (must be longer than 1.3653 s)")
end

% Check the channel number of the input data
if size(p, 2) > 2
    error("Error: Input signal comprises more than two channels")
else
    chansIn = size(p, 2);
    if chansIn == 2
        chans = ["Stereo left";
                 "Stereo right"];
    else
        chans = "Mono";
    end
end

%% Define constants

signalT = size(p, 1)/sampleRateIn;  % duration of input signal
sampleRate48k = 48e3;  % Signal sample rate prescribed to be 48kHz (to be used for resampling), Section 5.1.1 ECMA-418-2:2025 [r_s]

dz = 0.5;  % critical band overlap [deltaz]
halfBark = 0.5:dz:26.5;  % half-overlapping critical band rate scale [z]
nBands = length(halfBark);  % number of bands
bandCentreFreqs = shmBark2Hz(halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2025 [F(z)]

% Block and hop sizes Section 9.1.1 ECMA-418-2:2025
overlap = 0.75;  % block overlap proportion
blockSize = 65536;  % block size [s_b]
hopSize = (1 - overlap)*blockSize;  % hop size [s_h]

% Downsampled block and hop sizes Section 9.1.2 ECMA-418-2:2025
downSample = 32;  % downsampling factor
sampleRate1500 = sampleRate48k/downSample;  % [r~s]
blockSize1500 = blockSize/downSample;  % [s~b] = 2048
deltaF1500 = sampleRate1500/blockSize1500;  % DFT resolution Section 9.1.5 [Delta f]

% Determination of envelope analysis windows Section 9.1.3 ECMA-418-2:2025
num0WinStart = blockSize1500/32;  % start number of zeros [n_zb], initial value
num0WinEnd = num0WinStart;  % end number of zeros [n_ze], initial value
envWindow = zeros([blockSize1500, 1]);  % default envelope analysis window [w_Elz(ntilde)]
envWindow(num0WinStart + 1:end - num0WinEnd) = 1;
movMedLen = blockSize1500/64 + 1;  % moving median filter length (Section 9.1.3.2) = 33
quietZerosMin = blockSize1500*5/32;  % minimum duration of an interior quieter period (Section 9.1.3.4) [ntilde_zeros,min] = 320
winEndLim = blockSize1500/4;  % constraint on end index of envelope analysis window (Section 9.1.3.6), expressed as a MATLAB 1-based index (= ntilde_2 >= s~b/4 - 1 in the standard's 0-based indexing)
linRegIdxShift = blockSize1500*5/1024;  % shift on window indices for linear regression analysis (Section 9.1.3.6) = 10
hilbertMargin = blockSize1500/32;  % Hilbert-transform distortion margin added when updating window parameters (Section 9.1.3.5) = 64

% Section 9.1.4 - standardised epsilon (substituting for the standard's
% epsilon_0, the smallest positive double such that 1 + epsilon_0 > 1 -
% see the Note in shmHSA.m for why 1e-12 is used here without loss of
% accuracy, and for consistency with the repository's existing use of
% epsilon = 1e-12 throughout, e.g. in acousticSHMRoughness.m)
epsilon = 1e-12;

% term used to compensate for MATLAB 1-indexing
mlabIdx = 1;

% Section 9.1.5 - candidate modulation rates for the search for local
% minima of the HSA error function E_l,z((0,f_i)), f_i = 0.25*2^((i-2)/3)
% Hz, i = 1,...,16
modRateInitial = 0.25*2.^(((1:16) - 2)./3);
phiEmin = 0.15;  % Section 9.1.5 Equation 143 [Phi_Emin]

% Section 9.1.6 Equation 148 - fluctuation strength band-pass weighting
% is applied via shmFluctWeight.m using bandCentreFreqs(zBand) directly

% Section 9.1.7 - fine tuning (modified damped Newton method) parameters
newtonDx = 1e-5;  % finite-difference step [Delta x]
newtonMaxIt = 40;  % Equation 152 maximum number of iterations
% Equation 152 maximum step size [Hz]
%
% CORRECTION relative to the printed standard: Equation 152 as typeset in
% ECMA-418-2:2025 caps the Newton step at 2*10^-4 Hz. With the 1/4 damping
% factor and the 40-iteration limit of Equation 152, the fine tuning can
% then move the modulation rate by at most 40*0.25*2e-4 = 0.002 Hz, which
% (i) is far smaller than the resolution of the initial estimates (Delta f
% = 0.732 Hz for Equation 144, and the 1/3-octave grid of the f_min
% candidates), so the optimisation almost never converges and simply
% stops after 40 cap-limited steps, and (ii) makes the rejection test
% |f_c,1,opt - f_c,imax| > 1.25*Delta f = 0.92 Hz in Section 9.1.7
% unreachable. Diagnostics on real recordings confirmed the cap is active
% in 77-87 % of all blocks (fine-tuned rates clustering at exactly
% f_i +/- 0.002 Hz). With the cap raised to 2*10^-1 Hz the optimisation
% converges, the rejection test becomes meaningful (maximum travel 2 Hz),
% and the agreement of the time-dependent specific fluctuation strength
% with the reference implementation (HEAD acoustics ArtemiS v17) improves
% substantially on complex recordings (e.g. band-spectrum relative error
% reduced by about one third), with no effect on the calibration
% sinusoids (whose 4 Hz rate lies exactly on the candidate grid). The
% printed value is therefore treated as a typographical error in the
% exponent. As-written: newtonStepLim = 2e-4;
newtonStepLim = 2e-1;
newtonConvTol = 1e-7;  % Equation 152 convergence tolerance [Hz]
newtonRejectTol = 1.25*deltaF1500;  % Section 9.1.7 rejection tolerance

% Section 9.1.10 - threshold applied to A(l,z)
aThreshold = 5.2519;

% Output sample rate (Section 9.1.11 ECMA-418-2:2025) [r_s50]
sampleRate50 = 50;

% Calibration constant, Section 9.1.11 Equation 163 ECMA-418-2:2025 [c_F]
cal_F = 0.003840572;

%% Signal processing

% Input pre-processing
% --------------------
if sampleRateIn ~= sampleRate48k  % Resample signal
    [p_re, ~] = shmResample(p, sampleRateIn);
else  % don't resample
    p_re = p;
end

% Input signal samples
n_samples = size(p_re, 1);

% Section 5.1.2 ECMA-418-2:2025 Fade in weighting and zero-padding
% (only the start is zero-padded)
pn = shmPreProc(p_re, blockSize, hopSize, true, false);

% Apply outer & middle ear filter
% -------------------------------
%
% Section 5.1.3.2 ECMA-418-2:2025 Outer and middle/inner ear signal filtering
pn_om = shmOutMidEarFilter(pn, soundField);

n_steps = 53*4 + 6;  % approximate number of calculation steps, per channel

% Section 5.1.9 Table 3 ECMA-418-2:2025 - loudness threshold in quiet
LTQz = [0.3310, 0.1625, 0.1051, 0.0757, 0.0576, 0.0453, 0.0365, 0.0298,...
        0.0247, 0.0207, 0.0176, 0.0151, 0.0131, 0.0115, 0.0103, 0.0093,...
        0.0086, 0.0081, 0.0077, 0.0074, 0.0073, 0.0072, 0.0071, 0.0072,...
        0.0073, 0.0074, 0.0076, 0.0079, 0.0082, 0.0086, 0.0092, 0.0100,...
        0.0109, 0.0122, 0.0138, 0.0157, 0.0172, 0.0180, 0.0180, 0.0177,...
        0.0176, 0.0177, 0.0182, 0.0190, 0.0202, 0.0217, 0.0237, 0.0263,...
        0.0296, 0.0339, 0.0398, 0.0485, 0.0622];

% Loop through channels in file
% -----------------------------
if diagOn
    % cross-channel diagnostic log accumulators (see diagOn in the
    % function help and fluctuationSHM.diagLog below)
    diagLogAll = struct('chan', [], 'zBand', [], 'lBlock', [], 'status', {{}},...
                         'Mc', [], 'KL', [], 'nUnknown', [], 'rcondA', [],...
                         'usedPinv', [], 'fcOpt', [], 'nHarm', [], 'Ahat', [],...
                         'powerSum', [], 'N_HSA', [], 'A_lz', [], 'belowThreshold', [],...
                         'nzb', [], 'nze', []);
end

for chan = chansIn:-1:1

    if waitBar
        w = waitbar(0, "Initialising...");
        i_step = 1;

        waitbar(i_step/n_steps, w, 'Applying auditory filters...');
        i_step = i_step + 1;
    end % end of if branch for waitBar

    % Apply auditory filter bank
    % --------------------------
    % Filter equalised signal using 53 1/2-overlapping Bark filters
    % according to Section 5.1.4.2 ECMA-418-2:2025
    pn_omz = shmAuditoryFiltBank(pn_om(:, chan), false);

    % Note: At this stage, typical computer RAM limits impose a need to loop
    % through the critical bands rather than continue with a parallelised
    % approach, until later downsampling is applied
    for zBand = nBands:-1:1
        % Segmentation into blocks
        % ------------------------
        if waitBar
            waitbar(i_step/n_steps, w, strcat("Calculating signal envelopes in 53 bands, ",...
                           num2str(zBand), " to go..."));...
            i_step = i_step + 1;
        end % end of if branch for waitBar

        % Section 5.1.5 ECMA-418-2:2025
        i_start = 1;
        [pn_lz, iBlocksOut] = shmSignalSegment(pn_omz(:, zBand), 1,...
                                               blockSize, overlap,...
                                               i_start, true);

        % Envelope calculation and downsampling
        % --------------------------------------
        % Section 9.1.2 ECMA-418-2:2025
        % magnitude of Hilbert transform with downsample - Equation 119
        % [p(ntilde)_E,l,z]
        envelopes(:, :, zBand) = downsample(abs(hilbert(pn_lz)), downSample, 0);

    end  % end of for loop for obtaining low frequency signal envelopes

    % Note: With downsampled envelope signals, a parallelised approach can
    % continue for the window/spectrum determination

    nBlocks = size(envelopes, 2);

    % Determination of envelope analysis windows
    % -------------------------------------------
    % Section 9.1.3 ECMA-418-2:2025

    % Section 9.1.3.1 default window parameters, broadcast over blocks/bands
    num0WinStartMat = repmat(num0WinStart, nBlocks, nBands);
    num0WinEndMat = repmat(num0WinEnd, nBlocks, nBands);
    envWindowMat = repmat(envWindow, 1, nBlocks, nBands);

    % Section 9.1.3.2 Envelope smoothing and first weighting
    % moving median filter of length movMedLen, rounded to 8 decimal
    % places, windowed by the default window, then find the windowed max
    envMedWeight = envWindow.*round(movmedian(envelopes, movMedLen, 1), 8);  % [pBar(nTilde)_E,l,z] with initial window applied
    envMedWghtMax = max(envMedWeight, [], 1);  % [p_Emax,l,z]
    quietBlockInit = envMedWghtMax <= 5e-6;  % entire-block quieter period flag
    envWindowMat(:, quietBlockInit) = 0;

    if any(envWindowMat, 'all')
        % Section 9.1.3.3 Further quieter period detection
        quietThreshold = round(0.01*envMedWghtMax, 8);  % [p_Ethr,l,z]
        quietPeriods = envMedWeight < quietThreshold;

        % get the first and last indices of non-quiet samples (0-based
        % equivalents of n_zb and s~b - 1 - n_ze respectively)
        [iRows, jCols, kPages] = ind2sub([blockSize1500, nBlocks, nBands], find(~quietPeriods));
        totalBlocks = nBlocks*nBands;
        blockBandIdx = sub2ind([nBlocks, nBands], jCols, kPages);
        startIdx = reshape(accumarray(blockBandIdx, iRows, [totalBlocks, 1], @min, NaN), [nBlocks, nBands]);
        endIdx = reshape(accumarray(blockBandIdx, iRows, [totalBlocks, 1], @max, NaN), [nBlocks, nBands]);

        % update number of zeros (Section 9.1.3.3) - startIdx/endIdx are
        % 1-based MATLAB indices of the first/last non-quiet sample, so
        % (startIdx - 1) and (blockSize1500 - endIdx) are the 0-based
        % zero counts n_zb and n_ze respectively
        num0WinStartMatNew = max(num0WinStartMat, startIdx - 1, 'omitmissing');
        num0WinEndMatNew = max(num0WinEndMat, blockSize1500 - endIdx, 'omitmissing');
        startIdxNew = num0WinStartMatNew + 1;
        endIdxNew = blockSize1500 - num0WinEndMatNew;

        % Section 9.1.3.4-9.1.3.6 - search for the longest interior
        % quieter period, update the window accordingly, and validate
        % (these steps are adaptive per block/band and are therefore
        % looped, consistent with the equivalent adaptive steps in
        % acousticSHMRoughness.m)
        num0WinStartMatFinal = num0WinStartMatNew;
        num0WinEndMatFinal = num0WinEndMatNew;
        startIdxNewFinal = startIdxNew;
        endIdxNewFinal = endIdxNew;
        linRegStdChk = true(nBlocks, nBands);

        for zBand = 1:nBands
            for nBlock = 1:nBlocks
                if quietBlockInit(1, nBlock, zBand)
                    continue  % already flagged as an entire-block quieter period
                end

                % assign updated interval
                quietPeriodsNew = quietPeriods(startIdxNew(nBlock, zBand):endIdxNew(nBlock, zBand),...
                                               nBlock, zBand);

                % pad to identify start and end of quieter periods within
                % the updated interval
                quietPeriodsPad = [0; quietPeriodsNew; 0];

                quietPeriodStartsNew = find(diff(quietPeriodsPad(1:end - 1), 1) == 1);
                % Note: the subtraction of one is to identify the last
                % index of each quieter period, which is one before the
                % difference of -1 occurs
                quietPeriodEndsNew = find(diff(quietPeriodsPad(2:end), 1) == -1) - 1;

                quietPeriodLengths = quietPeriodEndsNew - quietPeriodStartsNew + 1;

                % Section 9.1.3.4 - minimum length criterion mask (a
                % quieter period duration must be STRICTLY greater than
                % quietZerosMin, per the standard's ">" in this clause)
                mask = quietPeriodLengths > quietZerosMin;
                if any(mask)
                    quietPeriodStartsNew = quietPeriodStartsNew(mask);
                    quietPeriodEndsNew = quietPeriodEndsNew(mask);
                    quietPeriodLengths = quietPeriodLengths(mask);

                    [~, qpLongestIdx] = max(quietPeriodLengths);

                    % get the start and end indices for the longest
                    % quieter period, adjusting for the updated start
                    % index [n~_qpmb,l,z], [n~_qpme,l,z] (1-based)
                    quietPeriodStartsMax = quietPeriodStartsNew(qpLongestIdx) + startIdxNew(nBlock, zBand) - 1;
                    quietPeriodEndsMax = quietPeriodEndsNew(qpLongestIdx) + startIdxNew(nBlock, zBand) - 1;

                    % Section 9.1.3.5 - update of the analysis window
                    % parameters. The standard's condition (using 0-based
                    % indices) is:
                    %   (n~qpmb - (nzb + 64)) > ((s~b - 1 - nze - 64) - n~qpme)
                    % which simplifies (the "+/-64" terms cancel) to a
                    % direct comparison of the left-hand and right-hand
                    % candidate active-window lengths. Working here with
                    % 1-based indices (quietPeriodStartsMax = n~qpmb + 1,
                    % quietPeriodEndsMax = n~qpme + 1), the equivalent
                    % comparison is:
                    leftLen = quietPeriodStartsMax - num0WinStartMatNew(nBlock, zBand) - 1;
                    rightLen = (blockSize1500 - num0WinEndMatNew(nBlock, zBand)) - quietPeriodEndsMax;

                    if leftLen > rightLen
                        % keep the left (longer) part: update n_ze only
                        num0WinEndMatFinal(nBlock, zBand) = blockSize1500 - quietPeriodStartsMax + hilbertMargin;
                    else
                        % keep the right (longer) part: update n_zb only
                        num0WinStartMatFinal(nBlock, zBand) = quietPeriodEndsMax - 1 + hilbertMargin;
                    end  % end of if branch to determine which window end is modified
                end  % end of if branch for minimum quieter period length

                % assign zeros to window start and end
                startIdxNewFinal(nBlock, zBand) = num0WinStartMatFinal(nBlock, zBand) + 1;
                endIdxNewFinal(nBlock, zBand) = blockSize1500 - num0WinEndMatFinal(nBlock, zBand);

                envWindowMat(1:num0WinStartMatFinal(nBlock, zBand), nBlock, zBand) = 0;
                envWindowMat(endIdxNewFinal(nBlock, zBand) + 1:end, nBlock, zBand) = 0;

                % Section 9.1.3.6 window interval validity checks
                % relative standard deviation of a linear regression fit
                if endIdxNewFinal(nBlock, zBand) - startIdxNewFinal(nBlock, zBand) + 1 >= quietZerosMin...
                        && endIdxNewFinal(nBlock, zBand) >= winEndLim
                    idxRange = (startIdxNewFinal(nBlock, zBand) + linRegIdxShift):(endIdxNewFinal(nBlock, zBand) - linRegIdxShift);
                    linReg = polyfit(idxRange, envMedWeight(idxRange, nBlock, zBand), 1);
                    pred = polyval(linReg, idxRange);
                    resid = pred.' - envMedWeight(idxRange, nBlock, zBand);
                    stdResid = std(resid);
                    meanVal = mean(envMedWeight(idxRange, nBlock, zBand));
                    linRegStdChk(nBlock, zBand) = meanVal ~= 0 && (stdResid/abs(meanVal)) >= 0.1/100;
                else
                    linRegStdChk(nBlock, zBand) = false;
                end
            end  % end of for loop over blocks for quieter periods
        end  % end of for loop over bands for quieter periods

        % Section 9.1.3.6 window interval validity checks
        % length of window
        winLengthChk = (endIdxNewFinal - startIdxNewFinal + 1) >= quietZerosMin;

        % window end is far enough into block
        winEndChk = endIdxNewFinal >= winEndLim;

        % all checks combined (reshape removes the leading singleton
        % dimension of quietBlockInit, giving an [nBlocks x nBands] array
        % that matches the shape of the other validity-check matrices;
        % reshape is used rather than squeeze so that this remains
        % correct even in the edge case nBlocks = 1)
        quietBlock = reshape(quietBlockInit, [nBlocks, nBands]) | ~(linRegStdChk & winLengthChk & winEndChk);

        num0WinStartMatFinal(quietBlock) = blockSize1500/2;
        num0WinEndMatFinal(quietBlock) = blockSize1500/2;

        envWindowMat(:, quietBlock) = 0;
    else
        quietBlock = true(nBlocks, nBands);
        num0WinStartMatFinal = repmat(blockSize1500/2, nBlocks, nBands);
        num0WinEndMatFinal = repmat(blockSize1500/2, nBlocks, nBands);
    end  % end of if branch for further quieter period analysis

    % High-resolution Spectral Analysis input spectra
    % -------------------------------------------------
    % Section 9.1.4 Equations 121-122 [P_E,l,z(k)], [Phi_E,l,z(k)]
    envSpectra = fft(envelopes.*envWindowMat, blockSize1500, 1);
    envMagSqSpectra = abs(envSpectra).^2;

    % Section 9.1.5 Equation 143 threshold [max(0.001*Phi_E,l,z(0), Phi_Emin)]
    modSpecCriterion = max(0.001*reshape(envMagSqSpectra(1, :, :), [nBlocks, nBands]), phiEmin);

    %% High-resolution Spectral Analysis pipeline (Sections 9.1.4-9.1.10)
    %
    % This part of the calculation is inherently adaptive (a
    % variable-dimension linear system is solved per block and per band,
    % with an iterative fine-tuning optimisation and harmonic-order
    % search), so - consistent with the equivalent adaptive steps in
    % acousticSHMRoughness.m (Sections 7.1.5.1 and 7.1.5.3) - it is
    % implemented with nested loops over bands and blocks rather than
    % parallelised across the whole array.

    AhatMat = zeros(nBlocks, nBands);        % [Ahat(l,z)], Equation 157
    powerSumMat = zeros(nBlocks, nBands);    % [phat_0^2 + 2*sum(A_i)], Equation 159 denominator (raw, unweighted power)
    N_HSA_Mat = zeros(nBlocks, nBands);      % [N'_HSA(l,z)], Equation 161
    fundRateMat = zeros(nBlocks, nBands);    % [f_1(l,z)], Equation 156

    % Diagnostic log preallocation (diagOn only - see diagOn in the
    % function help). Upper bound nBlocks*nBands rows; unused rows are
    % trimmed to diagIdx before being appended to the cross-channel log.
    if diagOn
        diagChan = zeros(nBlocks*nBands, 1);
        diagZBand = zeros(nBlocks*nBands, 1);
        diagLBlock = zeros(nBlocks*nBands, 1);
        diagStatus = cell(nBlocks*nBands, 1);
        diagMc = nan(nBlocks*nBands, 1);
        diagKL = nan(nBlocks*nBands, 1);
        diagNUnknown = nan(nBlocks*nBands, 1);
        diagRcondA = nan(nBlocks*nBands, 1);
        diagUsedPinv = nan(nBlocks*nBands, 1);
        diagFcOpt = nan(nBlocks*nBands, 1);
        diagNHarm = nan(nBlocks*nBands, 1);
        diagAhat = nan(nBlocks*nBands, 1);
        diagPowerSum = nan(nBlocks*nBands, 1);
        diagNHSA = nan(nBlocks*nBands, 1);
        diagAlz = nan(nBlocks*nBands, 1);
        diagBelowThresh = nan(nBlocks*nBands, 1);
        diagNzb = nan(nBlocks*nBands, 1);
        diagNze = nan(nBlocks*nBands, 1);
        diagIdx = 0;
    end

    for zBand = nBands:-1:1
        if waitBar
            waitbar(i_step/n_steps, w, strcat("Running HSA in 53 bands, ",...
                    num2str(zBand), " to go..."));...
            i_step = i_step + 1;
        end % end of if branch for waitBar

        nzbBand = num0WinStartMatFinal(:, zBand);
        nzeBand = num0WinEndMatFinal(:, zBand);

        for lBlock = 1:nBlocks
            if quietBlock(lBlock, zBand)
                continue  % A(l,z) remains zero for this block (quieter period)
            end

            nzb = nzbBand(lBlock);
            nze = nzeBand(lBlock);
            spectrumE = envSpectra(:, lBlock, zBand);
            spectrumPhi = envMagSqSpectra(1:49, lBlock, zBand);  % k = 0,...,48

            % Section 9.1.5 Stage 1 - local maxima of Phi_E,l,z(k) for
            % k = 1,...,47 (k = 0 and k = 48 cannot be local maxima given
            % the bounded search range, consistent with the need for
            % both neighbours k-1 and k+1 in Equation 144)
            [phiPks, kLocs] = findpeaks(spectrumPhi);
            keepMask = phiPks >= modSpecCriterion(lBlock, zBand);
            kLocs = kLocs(keepMask);  % 1-based MATLAB index; k0 = kLocs - mlabIdx
            % Section 9.1.5 - number of local maxima cannot exceed 24;
            % this is guaranteed by construction (at most floor(47/2)
            % alternating local maxima are possible over 47 interior
            % points) and is not separately enforced here.

            fpCandidates = zeros(size(kLocs));
            for iPk = 1:numel(kLocs)
                k0 = kLocs(iPk) - mlabIdx;  % 0-based peak index
                jIdx = (k0 - 1):(k0 + 1);  % 0-based neighbour indices
                phiNeighbours = spectrumPhi(jIdx + mlabIdx);
                % Section 9.1.5 Equation 144 [f_p,i(l,z)]
                %
                % CORRECTION relative to the printed standard: Equation 144
                % as typeset in ECMA-418-2:2025 includes a "- 1" term inside
                % the outer brackets, i.e. fp,i = (centroid - 1)*Delta_f.
                % This was verified to be present in the actual typeset
                % page image (not a text-extraction artefact), but applying
                % it produces a systematic bias of very close to one full
                % DFT bin (Delta_f) below the true frequency, confirmed
                % numerically across multiple independent single-tone test
                % cases (errors of -0.73 to -0.85 Hz, i.e. essentially
                % -Delta_f, versus +0.01 to -0.12 Hz - a small fraction of
                % one bin, consistent with a standard power-weighted
                % three-point centroid interpolator - with the "- 1" term
                % removed). The weighted centroid of the bin indices
                % themselves (without an additional offset) is the
                % conventional and correct form of this estimator, matching
                % the analogous (unbiased) refinement step used by
                % Equations 73-76 for roughness. Removing this term also
                % resolved a >90% amplitude recovery error for a secondary
                % (non-dominant) spectral line in a two-tone synthetic test,
                % which the biased estimate could push far enough from the
                % true frequency to substantially corrupt that line's HSA
                % fit; the dominant line was less affected only because
                % Section 9.1.5's own case I/II duplicate-selection logic
                % happened to discard the biased estimate in favour of
                % f_min in that specific (single-line) test case.
                % As-written: fpCandidates(iPk) = (sum(jIdx(:).*phiNeighbours)/sum(phiNeighbours) - 1)*deltaF1500;
                fpCandidates(iPk) = (sum(jIdx(:).*phiNeighbours)/sum(phiNeighbours))*deltaF1500;
            end

            % Section 9.1.5 Stage 2 - local minima of the HSA error
            % function E_l,z((0,f_i)) over the 16 log-spaced candidates.
            % Only interior points (i = 2,...,15) can be local minima;
            % this is what the standard's unqualified "local minima"
            % means for a bounded, ordered set of candidates.
            errHSA = zeros(1, 16);
            for iCand = 1:16
                [~, errHSA(iCand)] = shmHSA(modRateInitial(iCand), spectrumE,...
                                            blockSize1500, sampleRate1500, nzb, nze, epsilon);
            end
            isLocalMin = false(1, 16);
            isLocalMin(2:15) = errHSA(2:15) < errHSA(1:14) & errHSA(2:15) < errHSA(3:16);

            if ~any(isLocalMin)
                if isempty(fpCandidates)
                    % Section 9.1.5 - no local maximum and no local
                    % minimum: no modulation in this block
                    if diagOn
                        diagIdx = diagIdx + 1;
                        diagChan(diagIdx) = chan; diagZBand(diagIdx) = zBand; diagLBlock(diagIdx) = lBlock;
                        diagStatus{diagIdx} = 'no_modulation';
                        diagNzb(diagIdx) = nzb; diagNze(diagIdx) = nze;
                    end
                    continue
                end
                % Section 9.1.5 - no local minimum: use all local maxima
                fcFinal = sort(fpCandidates(:).');
                [pHatAll, ~, diagBig] = shmHSA(fcFinal, spectrumE, blockSize1500,...
                                     sampleRate1500, nzb, nze, epsilon);
            else
                errCandidates = errHSA(isLocalMin);
                freqCandidates = modRateInitial(isLocalMin);
                [~, iBestMin] = min(errCandidates);
                fMin = freqCandidates(iBestMin);

                % Section 9.1.5 Equation 145 - duplicate detection
                idDup = abs(fMin - fpCandidates) < 1.25*deltaF1500;

                if isempty(fpCandidates)
                    fcFinal = fMin;
                    [pHatAll, ~, diagBig] = shmHSA(fcFinal, spectrumE, blockSize1500,...
                                         sampleRate1500, nzb, nze, epsilon);
                elseif ~any(idDup)
                    fcFinal = sort([fpCandidates(:).', fMin]);
                    [pHatAll, ~, diagBig] = shmHSA(fcFinal, spectrumE, blockSize1500,...
                                         sampleRate1500, nzb, nze, epsilon);
                else
                    % Case I: fMin plus all local maxima except duplicates
                    fcCaseI = sort([fpCandidates(~idDup).', fMin]);
                    [pHatCaseI, errCaseI, diagCaseI] = shmHSA(fcCaseI, spectrumE, blockSize1500,...
                                                   sampleRate1500, nzb, nze, epsilon);

                    % Case II: all local maxima only
                    fcCaseII = sort(fpCandidates(:).');
                    [pHatCaseII, errCaseII, diagCaseII] = shmHSA(fcCaseII, spectrumE, blockSize1500,...
                                                     sampleRate1500, nzb, nze, epsilon);

                    if errCaseI <= errCaseII
                        fcFinal = fcCaseI;
                        pHatAll = pHatCaseI;
                        diagBig = diagCaseI;
                    else
                        fcFinal = fcCaseII;
                        pHatAll = pHatCaseII;
                        diagBig = diagCaseII;
                    end
                end
            end

            if isempty(fcFinal)
                if diagOn
                    diagIdx = diagIdx + 1;
                    diagChan(diagIdx) = chan; diagZBand(diagIdx) = zBand; diagLBlock(diagIdx) = lBlock;
                    diagStatus{diagIdx} = 'no_modulation';
                    diagMc(diagIdx) = diagBig.Mc; diagKL(diagIdx) = diagBig.KL;
                    diagNUnknown(diagIdx) = diagBig.nUnknown; diagRcondA(diagIdx) = diagBig.rcondA;
                    diagUsedPinv(diagIdx) = diagBig.usedPinv;
                end
                continue
            end

            % Section 9.1.5 Equation 146 - amplitude threshold (raw,
            % unweighted power A_i(l,z) = |P_HSA,i|^2; strictly greater
            % than, per the standard)
            % Note: shmHSA.m returns TWO-SIDED spectral line amplitudes
            % (half the cosine amplitude of each envelope component), so
            % that phat_0^2 + 2*sum(A_i) in Equations 159-160 is the
            % mean-square power of the harmonic complex - see the Note in
            % shmHSA.m regarding Equation 123 and footnote 46.
            aRaw = abs(pHatAll(2:end)).^2;
            keepFinal = aRaw > 0.05*max(aRaw);
            if ~any(keepFinal)
                if diagOn
                    diagIdx = diagIdx + 1;
                    diagChan(diagIdx) = chan; diagZBand(diagIdx) = zBand; diagLBlock(diagIdx) = lBlock;
                    diagStatus{diagIdx} = 'no_survivors';
                    diagMc(diagIdx) = diagBig.Mc; diagKL(diagIdx) = diagBig.KL;
                    diagNUnknown(diagIdx) = diagBig.nUnknown; diagRcondA(diagIdx) = diagBig.rcondA;
                    diagUsedPinv(diagIdx) = diagBig.usedPinv;
                    diagNzb(diagIdx) = nzb; diagNze(diagIdx) = nze;
                end
                continue
            end
            fcSurvive = fcFinal(keepFinal);
            aRawSurvive = aRaw(keepFinal);

            % Section 9.1.6 Equations 147-148 - weighted power spectrum
            % and dominant candidate [Atilde_i(l,z)], [i_max]
            wlh = shmFluctWeight(fcSurvive, bandCentreFreqs(zBand));
            aTildeSurvive = aRawSurvive(:).'.*wlh(:).';
            [~, iMax] = max(aTildeSurvive);

            % Section 9.1.7 Equations 149-152 - fine tuning of the
            % dominant modulation rate (modified damped Newton method)
            x0 = fcSurvive(iMax);
            xk = x0;
            for kIter = 1:newtonMaxIt
                [~, eMid] = shmHSA(xk, spectrumE, blockSize1500, sampleRate1500, nzb, nze, epsilon);
                [~, ePlus] = shmHSA(xk + newtonDx, spectrumE, blockSize1500, sampleRate1500, nzb, nze, epsilon);
                [~, eMinus] = shmHSA(xk - newtonDx, spectrumE, blockSize1500, sampleRate1500, nzb, nze, epsilon);

                dE = (ePlus - eMinus)/(2*newtonDx);       % Equation 149
                d2E = (ePlus - 2*eMid + eMinus)/newtonDx^2;  % Equation 150

                deltaX = 0.25*sign(dE)*min(abs(dE)/(abs(d2E) + epsilon), newtonStepLim);  % Equation 152
                xk = xk - deltaX;  % Equation 151

                if abs(deltaX) <= newtonConvTol
                    break
                end
            end
            fcOpt = xk;

            if abs(fcOpt - x0) > newtonRejectTol
                fcOpt = x0;  % optimisation rejected, retain original estimate
            else
                % Section 9.1.7 - replace the modulation rate of the
                % maximum with the fine-tuned value and update the
                % corresponding spectral component of P_HSA and
                % Atilde_imax accordingly (constant part plus one
                % spectral line pair, as used by the optimisation)
                fcSurvive(iMax) = fcOpt;
                [pHatOpt, ~] = shmHSA(fcOpt, spectrumE, blockSize1500,...
                                      sampleRate1500, nzb, nze, epsilon);
                aRawSurvive(iMax) = abs(pHatOpt(2))^2;
                aTildeSurvive(iMax) = aRawSurvive(iMax)*shmFluctWeight(fcOpt, bandCentreFreqs(zBand));
            end

            if fcOpt < 0.125
                if diagOn
                    diagIdx = diagIdx + 1;
                    diagChan(diagIdx) = chan; diagZBand(diagIdx) = zBand; diagLBlock(diagIdx) = lBlock;
                    diagStatus{diagIdx} = 'discarded_low_freq';
                    diagMc(diagIdx) = diagBig.Mc; diagKL(diagIdx) = diagBig.KL;
                    diagNUnknown(diagIdx) = diagBig.nUnknown; diagRcondA(diagIdx) = diagBig.rcondA;
                    diagUsedPinv(diagIdx) = diagBig.usedPinv; diagFcOpt(diagIdx) = fcOpt;
                    diagNzb(diagIdx) = nzb; diagNze(diagIdx) = nze;
                end
                continue  % Section 9.1.7 - modulation discarded
            end

            % Section 9.1.8 Equations 153-156 - harmonic analysis
            % Test assumed orders o = 1, 2, 3 of fcOpt
            bestEnergy = -Inf;
            bestIset = [];
            bestRatios = [];
            for order = 1:3
                fc1o = fcOpt/order;
                ratios = round(fcSurvive/fc1o);  % Equation 153
                ratios(ratios > 5) = 0;  % ratios greater than 5 excluded
                validRatio = ratios > 0;
                tolCheck = false(size(fcSurvive));
                tolCheck(validRatio) = abs(fcSurvive(validRatio)./(ratios(validRatio)*fc1o) - 1) < 0.04;  % Equation 154
                % Section 9.1.8 - a "harmonic complex with fundamental
                % modulation rate f_c,1,o" (Equation 154) is only taken to
                % exist if one of the components is itself at that
                % fundamental (integer ratio R = 1 within the 4 %
                % tolerance). For o = 1 this is always satisfied by
                % f_c,1,opt; for o = 2, 3 it requires a component near
                % f_c,1,opt/o. INTERPRETATION NOTE: the printed text does
                % not state this explicitly, but it mirrors the roughness
                % procedure of Section 7.1.5.3 (Equations 88-91), where
                % every candidate fundamental is itself one of the
                % detected components, and without it the o = 2, 3 index
                % sets (which admit all half- and third-integer multiples
                % of f_c,1,opt) almost always accumulate more energy than
                % the o = 1 set purely by admitting more members. On complex
                % recordings this reading, combined with the Equation 152
                % step-cap correction above, reduced the relative error of
                % the time-dependent specific fluctuation strength against
                % the reference implementation (ArtemiS v17) by about one
                % third to one half, whereas testing only o = 1, dropping the
                % harmonic complex, or dropping w_bw all made agreement
                % worse.
                if ~any(tolCheck) || ~any(tolCheck & ratios == 1)
                    continue
                end
                energyOrder = sum(aTildeSurvive(tolCheck));  % Equation 155
                if energyOrder > bestEnergy
                    bestEnergy = energyOrder;
                    bestIset = tolCheck;
                    bestRatios = ratios;
                    bestOrder = order;
                end
            end

            if isempty(bestIset)
                if diagOn
                    diagIdx = diagIdx + 1;
                    diagChan(diagIdx) = chan; diagZBand(diagIdx) = zBand; diagLBlock(diagIdx) = lBlock;
                    diagStatus{diagIdx} = 'no_harmonic_group';
                    diagMc(diagIdx) = diagBig.Mc; diagKL(diagIdx) = diagBig.KL;
                    diagNUnknown(diagIdx) = diagBig.nUnknown; diagRcondA(diagIdx) = diagBig.rcondA;
                    diagUsedPinv(diagIdx) = diagBig.usedPinv; diagFcOpt(diagIdx) = fcOpt;
                    diagNzb(diagIdx) = nzb; diagNze(diagIdx) = nze;
                end
                continue
            end

            % Section 9.1.8 - correct the modulation rates of the
            % retained components to exact integer multiples of the
            % fundamental (the component already equal to fcOpt is left
            % numerically unchanged by this correction, since its own
            % ratio is exactly bestOrder by construction)
            fc1Fund = fcOpt/bestOrder;  % Equation 156 [f_1(l,z)]
            fcHarmCorrected = bestRatios(bestIset)*fc1Fund;
            nHarm = numel(fcHarmCorrected);

            % Section 9.1.8 - re-run the HSA individually (constant part
            % plus one spectral line) for each retained harmonic, and
            % average the resulting constant-part estimates
            p0Estimates = zeros(nHarm, 1);
            aRawHarm = zeros(nHarm, 1);
            for iHarm = 1:nHarm
                [pHatHarm, ~] = shmHSA(fcHarmCorrected(iHarm), spectrumE, blockSize1500,...
                                       sampleRate1500, nzb, nze, epsilon);
                p0Estimates(iHarm) = real(pHatHarm(1));
                aRawHarm(iHarm) = abs(pHatHarm(2))^2;
            end
            wlhHarm = shmFluctWeight(fcHarmCorrected, bandCentreFreqs(zBand));
            aTildeHarm = aRawHarm(:).*wlhHarm(:);
            p0Final = mean(p0Estimates);

            % Section 9.1.9 Equations 157-158 - weighting the sum of the
            % harmonic complex
            sumATilde = sum(aTildeHarm);
            cog = sum(fcHarmCorrected(:).*aTildeHarm)/(sumATilde + epsilon);  % centre of gravity [Hz]
            wBw = 1 + 0.79577*abs(cog - fcOpt)^0.43461;  % Equation 158
            Ahat = wBw*sumATilde;  % Equation 157 [Ahat(l,z)]

            % Section 9.1.10 Equations 159-161 - power of the harmonic
            % complex (RAW, unweighted amplitudes - see the Note in the
            % HSA explanation regarding A_i vs Atilde_i) and HSA-based
            % loudness
            powerSum = p0Final^2 + 2*sum(aRawHarm);  % Equation 159/160 denominator & argument
            pRMS_HSA = sqrt(0.5*max(0, powerSum));
            N_tilde_HSA = shmLoudNonlin(pRMS_HSA);  % Equation 160
            if N_tilde_HSA >= LTQz(zBand)
                N_HSA = N_tilde_HSA - LTQz(zBand);  % Equation 161
            else
                N_HSA = 0;
            end

            AhatMat(lBlock, zBand) = Ahat;
            powerSumMat(lBlock, zBand) = powerSum;
            N_HSA_Mat(lBlock, zBand) = N_HSA;
            fundRateMat(lBlock, zBand) = fc1Fund;

            if diagOn
                diagIdx = diagIdx + 1;
                diagChan(diagIdx) = chan; diagZBand(diagIdx) = zBand; diagLBlock(diagIdx) = lBlock;
                diagStatus{diagIdx} = 'ok';
                diagMc(diagIdx) = diagBig.Mc; diagKL(diagIdx) = diagBig.KL;
                diagNUnknown(diagIdx) = diagBig.nUnknown; diagRcondA(diagIdx) = diagBig.rcondA;
                diagUsedPinv(diagIdx) = diagBig.usedPinv; diagFcOpt(diagIdx) = fcOpt;
                diagNHarm(diagIdx) = nHarm; diagAhat(diagIdx) = Ahat;
                diagPowerSum(diagIdx) = powerSum; diagNHSA(diagIdx) = N_HSA;
                diagNzb(diagIdx) = nzb; diagNze(diagIdx) = nze;
            end

        end  % end of for loop over blocks
    end  % end of for loop over bands

    % Section 9.1.10 Equation 159 - scaling with HSA-based loudness,
    % completed here (requires max_z(N'_HSA(l,z)) across all bands for
    % each block, hence deferred to this second pass)
    N_HSA_max = max(N_HSA_Mat, [], 2);  % [max_z(N'_HSA(l,z))], per block
    A_lz = AhatMat.*(N_HSA_Mat.^2)./(N_HSA_max + epsilon)./(powerSumMat + epsilon).*blockSize1500;

    % threshold: values of A(l,z) below 5.2519 are set to zero, and the
    % corresponding fundamental modulation rates are also zeroed
    belowThreshold = A_lz < aThreshold;
    A_lz(belowThreshold) = 0;
    fundRateMat(belowThreshold) = 0; %#ok<NASGU> % retained for completeness/diagnostics; not used further below

    % Diagnostic log: backfill the pre-threshold A(l,z) value and the
    % threshold outcome for every logged block/band (only meaningful for
    % status = 'ok' rows; NaN elsewhere since A_lz is only computed - and
    % only nonzero prior to thresholding - for blocks that reached that
    % point in the pipeline)
    if diagOn
        linIdx = sub2ind([nBlocks, nBands], diagLBlock(1:diagIdx), diagZBand(1:diagIdx));
        isOkRow = strcmp(diagStatus(1:diagIdx), 'ok');
        % recomputed from the per-block matrices rather than read from
        % A_lz directly, since A_lz has already been zeroed above for
        % below-threshold entries and the pre-threshold value is the
        % diagnostically useful one; only meaningful (non-NaN) for rows
        % that actually reached the point where these matrices are
        % populated, i.e. status = 'ok'
        alzRecomputed = AhatMat(linIdx).*(N_HSA_Mat(linIdx).^2)./(N_HSA_max(diagLBlock(1:diagIdx)) + epsilon)./(powerSumMat(linIdx) + epsilon).*blockSize1500;
        diagAlz(1:diagIdx) = alzRecomputed;
        diagAlz(~isOkRow) = NaN;
        diagBelowThresh(1:diagIdx) = double(belowThreshold(linIdx));
        diagBelowThresh(~isOkRow) = NaN;
    end

    % Time-dependent specific fluctuation strength
    % ---------------------------------------------
    % Section 9.1.11 ECMA-418-2:2025

    % interpolation to 50 Hz sampling rate
    % Section 9.1.11 Equation 162 [t(l)] (via iBlocksOut from
    % shmSignalSegment with endShrink = true, which already places the
    % final block at the true end of the signal, matching Equation 162's
    % l = l_last special case)
    l_50Last = floor(n_samples/sampleRate48k*sampleRate50) + 1;
    x = (iBlocksOut - 1)/sampleRate48k;
    xq = linspace(0, signalT, l_50Last);
    specFluctEst = zeros(l_50Last, nBands);
    for zBand = nBands:-1:1
        specFluctEst(:, zBand) = pchip(x, A_lz(:, zBand), xq);
    end  % end of for loop for interpolation
    specFluctEst(specFluctEst < 0) = 0;  % [F'_est(l_50,z)]

    % Section 9.1.11 Equations 166-167 [Ftilde'_est(l_50)], [Fbar'_est(l_50)]
    specFluctEstRMS = rms(specFluctEst, 2);
    specFluctEstAvg = mean(specFluctEst, 2);

    % Section 9.1.11 Equation 165 [Bhat(l_50)], smoothed with a moving
    % median filter of length 71 to give [B(l_50)]
    Bhatl50 = zeros(size(specFluctEstAvg));
    mask = specFluctEstAvg ~= 0;
    Bhatl50(mask) = specFluctEstRMS(mask)./specFluctEstAvg(mask);
    Bl50 = movmedian(Bhatl50, 71);

    % Section 9.1.11 Equation 164 [E(l_50)]
    El50 = 0.37106*(tanh(1.6407*(Bl50 - 2.5804)) + 1)*0.5 + 0.58449;

    % Section 9.1.11 Equation 163 [Fhat'(l_50,z)]
    specFluctEstTform = cal_F*(specFluctEst.^El50);

    % Section 9.1.11 Equations 168 [F'(l_50,z)] - single time constant
    % low-pass filter of order one (equal rise/fall time constants, so
    % shmRoughLowPass.m - which implements the structurally identical
    % Equations 109-110 for roughness - is reused directly)
    tau = 0.75;
    specFluctuation(:, :, chan) = shmRoughLowPass(specFluctEstTform, sampleRate50, tau, tau);

    if waitBar
        close(w)  % close waitbar
    end

    if diagOn
        diagLogAll.chan = [diagLogAll.chan; diagChan(1:diagIdx)];
        diagLogAll.zBand = [diagLogAll.zBand; diagZBand(1:diagIdx)];
        diagLogAll.lBlock = [diagLogAll.lBlock; diagLBlock(1:diagIdx)];
        diagLogAll.status = [diagLogAll.status; diagStatus(1:diagIdx)];
        diagLogAll.Mc = [diagLogAll.Mc; diagMc(1:diagIdx)];
        diagLogAll.KL = [diagLogAll.KL; diagKL(1:diagIdx)];
        diagLogAll.nUnknown = [diagLogAll.nUnknown; diagNUnknown(1:diagIdx)];
        diagLogAll.rcondA = [diagLogAll.rcondA; diagRcondA(1:diagIdx)];
        diagLogAll.usedPinv = [diagLogAll.usedPinv; diagUsedPinv(1:diagIdx)];
        diagLogAll.fcOpt = [diagLogAll.fcOpt; diagFcOpt(1:diagIdx)];
        diagLogAll.nHarm = [diagLogAll.nHarm; diagNHarm(1:diagIdx)];
        diagLogAll.Ahat = [diagLogAll.Ahat; diagAhat(1:diagIdx)];
        diagLogAll.powerSum = [diagLogAll.powerSum; diagPowerSum(1:diagIdx)];
        diagLogAll.N_HSA = [diagLogAll.N_HSA; diagNHSA(1:diagIdx)];
        diagLogAll.A_lz = [diagLogAll.A_lz; diagAlz(1:diagIdx)];
        diagLogAll.belowThreshold = [diagLogAll.belowThreshold; diagBelowThresh(1:diagIdx)];
        diagLogAll.nzb = [diagLogAll.nzb; diagNzb(1:diagIdx)];
        diagLogAll.nze = [diagLogAll.nze; diagNze(1:diagIdx)];
    end

    clearvars envelopes  % avoid stale data carried over between channels

end  % end of for loop over channels

% Binaural fluctuation strength
% Section 9.1.15 ECMA-418-2:2025 [F'_B(l_50,z)]
if chansIn == 2 && binaural
    % Equation 170
    specFluctuation(:, :, 3) = sqrt(sum(specFluctuation(:, :, 1:2).^2, 3)/2);
    chansOut = 3;  % set number of 'channels' to stereo plus single binaural
    chans = [chans;
             "Binaural"];
else
    chansOut = chansIn;  % assign number of output channels
end

% Section 9.1.12 ECMA-418-2:2025
% Time-averaged specific fluctuation strength [F'(z)], discarding
% 0 <= l50 <= 35
specFluctuationAvg = mean(specFluctuation(37:end, :, :), 1);

% Section 9.1.13 ECMA-418-2:2025
% Time-dependent fluctuation strength Equation 169 [F(l_50)]
% Discard singleton dimensions
if chansOut == 1
    fluctuationTDep = sum(specFluctuation.*dz, 2);
    specFluctuationAvg = transpose(specFluctuationAvg);
else
    fluctuationTDep = squeeze(sum(specFluctuation.*dz, 2));
    specFluctuationAvg = squeeze(specFluctuationAvg);
end

% Section 9.1.14 ECMA-418-2:2025
% Overall fluctuation strength [F], discarding 0 <= l50 <= 35
fluctuation90Pc = prctile(fluctuationTDep(37:end, :), 90, 1);

% time (s) corresponding with results output [t]
timeOut = transpose((0:(size(specFluctuation, 1) - 1))/sampleRate50);

%% Output plotting

if outPlot
    % Plot figures
    % ------------
    for chan = chansOut:-1:1
        % Plot results
        fig = figure;
        tiledlayout(fig, 2, 1);
        movegui(fig, 'center');
        ax1 = nexttile(1);
        surf(ax1, timeOut, bandCentreFreqs, permute(specFluctuation(:, :, chan),...
                                              [2, 1, 3]),...
             'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        ax1.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        ax1.YLim = [bandCentreFreqs(1), bandCentreFreqs(end)];
        ax1.CLim = [0, max(1e-6, ceil(max(specFluctuation(:, :, chan), [], 'all')*500)/500)];
        ax1.YTick = [63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3]; 
        ax1.YTickLabel = ["63", "125", "250", "500", "1k", "2k", "4k",...
                          "8k", "16k"];
        ax1.YScale = 'log';
        ax1.YLabel.String = 'Frequency, Hz';
        ax1.XLabel.String = 'Time, s';
        ax1.FontName =  'Arial';
        ax1.FontSize = 12;
        cmap_inferno = load('cmap_inferno.txt');
        colormap(cmap_inferno);
        h = colorbar;
        set(get(h,'label'),'string', {'Specific fluctuation strength,'; 'vacil_{HMS}/Bark_{HMS}'});
        chan_lab = chans(chan);

        % Create A-weighting filter
        weightFilt = weightingFilter('A-weighting', sampleRateIn);
        % Filter signal to determine A-weighted time-averaged level
        if chan == 3
            pA = weightFilt(p);
            LAeq2 = 20*log10(rms(pA, 1)/2e-5);
            % take the higher channel level as representative (PD ISO/TS
            % 12913-3:2019 Annex D)
            [LAeq, LR] = max(LAeq2);
            % if branch to identify which channel is higher
            if LR == 1
                whichEar = ' left ear';
            else
                whichEar = ' right ear';
            end  % end of if branch

            chan_lab = chan_lab + whichEar;

        else
            pA = weightFilt(p(:, chan));
            LAeq = 20*log10(rms(pA)/2e-5);
        end
        
        title(strcat(chan_lab,...
                     ' signal sound pressure level =', {' '},...
                     num2str(round(LAeq,1)), "dB {\itL}_{Aeq}"),...
                     'FontWeight', 'normal', 'FontName', 'Arial');

        ax2 = nexttile(2);
        plot(ax2, timeOut, fluctuation90Pc(1, chan)*ones(size(timeOut)), ':', 'color',...
             cmap_inferno(34, :), 'LineWidth', 1.5, 'DisplayName', "90th" + string(newline) + "percentile");
        hold on
        plot(ax2, timeOut, fluctuationTDep(:, chan), 'color', cmap_inferno(166, :),...
             'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold off
        ax2.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        if max(fluctuationTDep(:, chan)) > 0
            ax2.YLim = [0, 1.1*ceil(max(fluctuationTDep(:, chan))*10)/10];
        end
        ax2.XLabel.String = 'Time, s';
        ax2.YLabel.String = 'Fluctuation strength, vacil_{HMS}';
        ax2.XGrid = 'on';
        ax2.YGrid = 'on';
        ax2.GridAlpha = 0.075;
        ax2.GridLineStyle = '--';
        ax2.GridLineWidth = 0.25;
        ax2.FontName = 'Arial';
        ax2.FontSize = 12;
        lgd = legend('Location', 'eastoutside', 'FontSize', 8);
        lgd.Title.String = "Overall";
    end  % end of for loop for plotting over channels
end  % end of if branch for plotting

%% Output assignment

% Assign outputs to structure
if chansOut == 3
    fluctuationSHM.specFluctuation = specFluctuation(:, :, 1:2);
    fluctuationSHM.specFluctuationAvg = specFluctuationAvg(:, 1:2);
    fluctuationSHM.fluctuationTDep = fluctuationTDep(:, 1:2);
    fluctuationSHM.fluctuation90Pc = fluctuation90Pc(:, 1:2);
    fluctuationSHM.specFluctuationBin = specFluctuation(:, :, 3);
    fluctuationSHM.specFluctuationAvgBin = specFluctuationAvg(:, 3);
    fluctuationSHM.fluctuationTDepBin = fluctuationTDep(:, 3);
    fluctuationSHM.fluctuation90PcBin = fluctuation90Pc(:, 3);
    fluctuationSHM.bandCentreFreqs = bandCentreFreqs;
    fluctuationSHM.timeOut = timeOut;
    fluctuationSHM.soundField = soundField;
else
    fluctuationSHM.specFluctuation = specFluctuation;
    fluctuationSHM.specFluctuationAvg = specFluctuationAvg;
    fluctuationSHM.fluctuationTDep = fluctuationTDep;
    fluctuationSHM.fluctuation90Pc = fluctuation90Pc;
    fluctuationSHM.bandCentreFreqs = bandCentreFreqs;
    fluctuationSHM.timeOut = timeOut;
    fluctuationSHM.soundField = soundField;
end

if diagOn
    fluctuationSHM.diagLog = diagLogAll;
end

% end of function