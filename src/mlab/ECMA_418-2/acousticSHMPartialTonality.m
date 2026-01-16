function partTonalitySHM = acousticSHMPartialTonality(pTarget, pMasker, sampleRateIn, axisN, soundField, waitBar, outPlot, annoyWeight)
% partTonalitySHM = acousticSHMPartialTonality(pTarget, pMasker, sampleRateIn, axisN, soundField, waitBar, outPlot, annoyWeight)
%
% Returns partial tonality values and frequencies based on ECMA-418-2:2025
% (using the Sottek Hearing Model) for an input calibrated single mono or single stereo
% audio (sound pressure) time-series target signal, pTarget, against a
% corresponding audio time-series masking signal, pMasker.
%
% Inputs
% ------
% pTarget : vector or 2D matrix
%   Input target signal as single mono or stereo audio (sound
%   pressure) signals (size must match pMasker).
%
% pMasker : vector or 2D matrix
%   Input masker signal as single mono or stereo audio (sound
%   pressure) signals (size must match pTarget).
%
% sampleRateIn : integer
%   Sample rate (frequency) of the input signals.
%
% axisN : integer (1 or 2, default: 1)
%   Time axis along which to calculate the tonality.
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
%   Determines whether a progress bar displays during processing.
%
% outPlot : Boolean true/false (default: false)
%   Flag indicating whether to generate a figure from the output
%   (set outPlot to false for doing multi-file parallel calculations).
%
% annoyWeight : Boolean true/false (default: false)
%   Flag indicating whether to include tonal annoyance-weighted results in
%   the output. This weighting is not currently included in the ECMA-418-2
%   standard, but was initially proposed by Sottek & Becker (2019).
% 
% Returns
% -------
%
% partTonalitySHM : structure
%   Contains the output.
%
% partTonalitySHM contains the following outputs:
%
% specPartTonality : matrix
%   Time-dependent specific partial tonality for each critical band
%   arranged as [time, bands(, channels)]
%
% specPartTonalityFreqs : matrix
%   Time-dependent frequencies of the dominant tonal
%   components corresponding with each of the time-dependent specific
%   partial tonality values in each critical band
%   arranged as [time, bands(, channels)]
%
% specPartTonalityAvg : matrix
%   Time-averaged specific partial tonality for each critical band
%   arranged as [bands(, channels)]
%
% specPartTonalityAvgFreqs : matrix
%   Frequencies of the dominant tonal components
%   corresponding with each of the
%   time-averaged specific partial tonality values in each critical band
%   arranged as [bands(, channels)]
%
% specPartTonalLoudness : matrix
%   Time-dependent specific tonal loudness for each critical band
%   arranged as [time, bands(, channels)]
%
% specNoiseLoudness : matrix
%   Time-dependent specific noise loudness for each critical band
%   arranged as [time, bands(, channels)]
%
% partTonalityTDep : vector or matrix
%   Time-dependent overall partial tonality arranged as [time(, channels)]
%
% partTonalityTDepFreqs : vector or matrix
%   Time-dependent frequencies of the dominant tonal components
%   corresponding with the time-dependent overall partial tonality values
%   arranged as [time(, channels)]
%
% partTonalityAvg : number or vector
%   Time-averaged overall partial tonality arranged as [tonality(, channels)]
%
% bandCentreFreqs : vector
%   Centre frequencies corresponding with each critical band rate
%
% timeOut : vector
%   Time (seconds) corresponding with time-dependent outputs
%
% soundField : string
%   Identifies the soundfield type applied (= input argument)
%
% If annoyWeight=true, additional outputs are included in partial tonality:
%
% specPartTonalityAnnoy : matrix
%   Tonal annoyance-weighted, time-dependent specific partial tonality for each
%   critical band arranged as [time, bands(, channels)]
%
% specPartTonalityAvgAnnoy : matrix
%   Tonal annoyance-weighted, time-averaged specific partial tonality for each
%   critical band arranged as [bands(, channels)]
%
% partTonalityTDepAnnoy : vector or matrix
%   Tonal annoyance-weighted, time-dependent overall partial tonality
%   arranged as [time(, channels)]
%
% partTonalityAvgAnnoy : number or vector
%   Tonal annoyance-weighted, time-averaged overall partial tonality
%   arranged as [tonality(, channels)]
%
% If outPlot=true, a set of plots is returned illustrating the energy
% time-averaged A-weighted sound level, the time-dependent specific and
% overall partial tonality, with the latter also indicating the time-aggregated
% value. A set of plots is returned for each input channel.
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
% References
% ----------
% Sottek, R & Becker, J (2019). Tonal annoyance vs. tonal loudness and
% tonality. In: Proceedings of Inter-noise 2019, Madrid, Spain, 13-16 June
% 2019.
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk) &
%          Matt Torjussen (m.c.torjussen@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 16/01/2026
% Date last modified: 16/01/2026
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
% This code calls sub-component file 'cmap_plasma.txt'. The contents of
% the file includes a copy of data obtained from the repository 
% https://github.com/BIDS/colormap, and is CC0 1.0 licensed for modified
% use, see https://creativecommons.org/publicdomain/zero/1.0 for
% information.
%
% This code was developed from an original file 'SottekTonality.m' authored
% by Matt Torjussen (14/02/2022), based on implementing ECMA-418-2:2020.
% The original code has been reused and updated here with permission.
%
% Checked by:
% Date last checked:
%
%% Arguments validation
    arguments (Input)
        pTarget (:, :) double {mustBeReal}
        pMasker (:, :) double {mustBeReal}
        sampleRateIn (1, 1) double {mustBePositive, mustBeInteger}
        axisN (1, 1) {mustBeInteger, mustBeInRange(axisN, 1, 2)} = 1
        soundField (1, :) string {mustBeMember(soundField,...
                                               {'freeFrontal',...
                                                'diffuse', ...
                                                'noOuter', ...
                                                'noEar'})} = 'freeFrontal'
        waitBar {mustBeNumericOrLogical} = true
        outPlot {mustBeNumericOrLogical} = false
        annoyWeight {mustBeNumericOrLogical} = false
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Input checks
% Target and masker signal size check
if size(pTarget) ~= size(pMasker)
    error("Error: Input target and masker signals must match in size.")
end

% Orient input matrix
if axisN == 2
    pTarget = pTarget.';
    pMasker = pMasker.';
end

% Check the length of the input data (must be longer than 300 ms)
if size(pTarget, 1) <  300/1000*sampleRateIn
    error("Error: Input signal is too short along the specified axis to calculate partial loudness (must be longer than 300 ms)")
end

% Check the channel number of the input data
if size(pTarget, 2) > 2
    error('Error: Input signal comprises more than two channels')
else
    chansIn = size(pTarget, 2);
    if chansIn > 1
        chans = ["Stereo left";
                 "Stereo right"];
    else
        chans = "Mono";
    end
end

%% Define constants

sampleRate48k = 48e3;  % Signal sample rate prescribed to be 48kHz (to be used for resampling), Section 5.1.1 ECMA-418-2:2025 [r_s]
deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2025 [deltaf(f=0)]
c = 0.1618;  % Half-overlapping Bark band centre-frequency denominator constant defined in Section 5.1.4.1 ECMA-418-2:2025

dz = 0.5;  % critical band overlap [deltaz]
halfBark = 0.5:dz:26.5;  % half-overlapping critical band rate scale [z]
bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2025 [F(z)]
dfz = sqrt(deltaFreq0^2 + (c*bandCentreFreqs).^2);  % Section 5.1.4.1 Equation 10 ECMA-418-2:2025 [deltaf(z)]

% Block and hop sizes Section 6.2.2 Table 4 ECMA-418-2:2025
overlap = 0.75;  % block overlap proportion
% block sizes [s_b(z)]
blockSize = [8192*ones(1, 3), 4096*ones(1, 13), 2048*ones(1, 9), 1024*ones(1, 28)];
% hop sizes (section 5.1.2 footnote 3 ECMA 418-2:2022) [s_h(z)]
hopSize = (1 - overlap)*blockSize;

% Output sample rate based on hop sizes - Resampling to common time basis
% Section 6.2.6 ECMA-418-2:2025 [r_sd]
sampleRate1875 = sampleRate48k/min(hopSize);

% Number of bands that need averaging. Section 6.2.3 Table 5
% ECMA-418-2:2025 [NB]
NBandsAvg = [0, 1, 2*ones(1,14), ones(1,9), zeros(1,28);...
             1, 1, 2*ones(1,14), ones(1,9), zeros(1,28)];

% Critical band interpolation factors from Section 6.2.6 Table 6
% ECMA-418-2:2025 [i]
i_interp = blockSize/min(blockSize);

% Noise reduction constants from Section 6.2.7 Table 7 ECMA-418-2:2025
alpha = 20;
beta = 0.07;

% Sigmoid function factor parameters Section 6.2.7 Table 8 ECMA-418-2:2025
% [c(s_b(z))]
csz_b = [18.21*ones(1, 3), 12.14*ones(1, 13), 417.54*ones(1, 9),...
         962.68*ones(1, 28)]; 
% [d(s_b(z))]
dsz_b = [0.36*ones(1, 3), 0.36*ones(1, 13), 0.71*ones(1, 9),...
         0.69*ones(1, 28)]; 

% Scaling factor constants from Section 6.2.8 Table 9 ECMA-418-2:2025
A = 35;
B = 0.003;

cal_T = 2.8758615;  % calibration factor in Section 6.2.8 Equation 51 ECMA-418-2:2025 [c_T]
cal_Tx = 1/0.9999043734252;  % Adjustment to calibration factor (Footnote 22 ECMA-418-2:2025)

% standardised epsilon
epsilon = 1e-12;

%% Signal processing

% Input pre-processing
% --------------------
if sampleRateIn ~= sampleRate48k  % Resample signal
    [pTarget_re, ~] = shmResample(pTarget, sampleRateIn);
    [pMasker_re, ~] = shmResample(pMasker, sampleRateIn);
else  % don't resample
    pTarget_re = pTarget;
    pMasker_re = pMasker;
end

% memory cleanup
clear pTarget pMAsker

% Section 5.1.2 ECMA-418-2:2025 Fade in weighting and zero-padding
pTn = shmPreProc(pTarget_re, max(blockSize), max(hopSize));
pMn = shmPreProc(pMasker_re, max(blockSize), max(hopSize));

% Apply outer & middle ear filter
% -------------------------------
%
% Section 5.1.3.2 ECMA-418-2:2025 Outer and middle/inner ear signal filtering
pTn_om = shmOutMidEarFilter(pTn, soundField);
pMn_om = shmOutMidEarFilter(pMn, soundField);

n_steps = 115;  % approximate number of calculation steps

% Loop through channels in file
% -----------------------------
for chan = chansIn:-1:1

    % Apply auditory filter bank
    % --------------------------

    if waitBar
        w = waitbar(0, "Initialising...");
        i_step = 1;

        waitbar(i_step/n_steps, w, 'Applying auditory filters...');
        i_step = i_step + 1;

    end % end of if branch for waitBar

    % Filter equalised signal using 53 1/2-overlapping Bark filters
    % according to Section 5.1.4.2 ECMA-418-2:2025
    pTn_omz = shmAuditoryFiltBank(pTn_om(:, chan));
    pMn_omz = shmAuditoryFiltBank(pMn_om(:, chan));

    % Autocorrelation function analysis
    % ---------------------------------
    % Duplicate Banded Data for ACF
    % Averaging occurs over neighbouring bands, to do this the segmentation
    % needs to be duplicated for neigbouring bands. 'Dupe' has been added
    % to variables to indicate that the vectors/matrices have been modified
    % for duplicated neigbouring bands.

    pTn_omzDupe = [pTn_omz(:, 1:5), pTn_omz(:, 2:18), pTn_omz(:, 16:26),...
                   pTn_omz(:, 26:53)];

    pMn_omzDupe = [pMn_omz(:, 1:5), pMn_omz(:, 2:18), pMn_omz(:, 16:26),...
                   pMn_omz(:, 26:53)];

    clear pTn_omz pMn_omz

    blockSizeDupe = [8192*ones(1, 5), 4096*ones(1, 17), 2048*ones(1, 11),...
                     1024*ones(1, 28)];
    bandCentreFreqsDupe = [bandCentreFreqs(1:5),...
                           bandCentreFreqs(2:18),...
                           bandCentreFreqs(16:26),...
                           bandCentreFreqs(26:53)];

    % (duplicated) indices corresponding with the NB bands around each z band
    i_NBandsAvgDupe = [1, 1, 1, 6:18, 23:31, 34:61;
                       2, 3, 5, 10:22, 25:33, 34:61];

    for zBand = 61:-1:1
        if waitBar
            waitbar(((62 - zBand) + i_step)/n_steps, w, strcat("Applying ACF in 61 bands, ",...
                num2str(zBand), " to go..."));
        end % end of if branch for waitBar

        % Segmentation into blocks
        % ------------------------
        % Section 5.1.5 ECMA-418-2:2025
        i_start = blockSizeDupe(1) - blockSizeDupe(zBand) + 1;
        [pTn_lz, ~] = shmSignalSegment(pTn_omzDupe(:, zBand), 1,...
                                       blockSizeDupe(zBand), overlap, i_start);
        [pMn_lz, ~] = shmSignalSegment(pMn_omzDupe(:, zBand), 1,...
                                       blockSizeDupe(zBand), overlap, i_start);

        % Transformation into Loudness
        % ----------------------------
        % Sections 5.1.6 to 5.1.9 ECMA-418-2:2025 [N'_basis(z)]
        [pTn_rlz, bandPartBasisLoudness, ~]...
            = shmBasisPartialLoudness(pTn_lz, pMn_lz, bandCentreFreqsDupe(zBand));

        % Apply ACF
        % ACF implementation using DFT
        % Section 6.2.2 Equations 27 & 28 ECMA-418-2:2025
        % [phi_unscaled,l,z(m)]
        unscaledACF = ifft(abs(fft(pTn_rlz, 2*blockSizeDupe(zBand), 1)).^2,...
                           2*blockSizeDupe(zBand), 1);
        % Section 6.2.2 Equation 29 ECMA-418-2:2025 [phi_l,z(m)]
        denom = sqrt(cumsum(pTn_rlz.^2, 1, 'reverse').*flipud(cumsum(pTn_rlz.^2, 1)))...
                + epsilon;

        % note that the block length is used here, rather than the 2*s_b,
        % for compatability with the remaining code - beyond 0.75*s_b is
        % assigned (unused) zeros in the next line
        unbiasedNormACF = unscaledACF(1:blockSizeDupe(zBand), :)./denom;
        unbiasedNormACF((0.75*blockSizeDupe(zBand) + 1):blockSizeDupe(zBand), :) = 0;

        % Section 6.2.2 Equation 30 ECMA-418-2:2025 [phi_z'(m)
        unbiasedNormACFDupe{zBand} = bandPartBasisLoudness.*unbiasedNormACF;

    end

    if waitBar
        i_step = i_step + 62;  % increment calculation step for waitbar
    end

    % Average the ACF over nB bands - Section 6.2.3 ECMA-418-2:2025        
    for zBand = 53:-1:1  % Loop through 53 critical band filtered signals
        if waitBar
            waitbar(((54 - zBand) + i_step)/n_steps, w,...
                    strcat("Calculating sound quality in 53 bands, ",...
                           num2str(zBand), " to go..."));
        end % end of if branch for waitBar

        NBZ = NBandsAvg(1, zBand) + NBandsAvg(2, zBand) + 1; % Total number of bands to average over

        % Averaging of frequency bands
        meanScaledACF = mean(reshape(cell2mat(unbiasedNormACFDupe(i_NBandsAvgDupe(1, zBand):i_NBandsAvgDupe(2, zBand))),...
                                     blockSize(zBand), [], NBZ), 3);

        % Average the ACF over adjacent time blocks [phibar_z'(m)]
        if zBand <= 16 
            meanScaledACF(:, 2:end-1) = movmean(meanScaledACF, 3, 2, 'omitnan',...
                                                'EndPoints', 'discard');
        end

        % Application of ACF lag window Section 6.2.4 ECMA-418-2:2025
        tauz_start = max(0.5/dfz(zBand), 2e-3);  % Equation 31 ECMA-418-2:2025 [tau_start(z)]
        tauz_end = max(4/dfz(zBand), tauz_start + 1e-3);  % Equation 32 ECMA-418-2:2025 [tau_end(z)]
        % Equations 33 & 34 ECMA-418-2:2025
        mz_start = ceil(tauz_start*sampleRate48k);  % Starting lag window index [m_start(z)]
        mz_end = floor(tauz_end*sampleRate48k);  % Ending lag window index [m_end(z)]
        M = mz_end - mz_start + 1;
        % Equation 35 ECMA-418-2:2025
        % lag-windowed, detrended ACF [phi'_z,tau(m)]
        lagWindowACF = zeros(size(meanScaledACF));
        lagWindowACF(mz_start:mz_end, :) = meanScaledACF(mz_start:mz_end, :)...
                                           - mean(meanScaledACF(mz_start:mz_end, :));

        % Estimation of partial tonal loudness
        % ------------------------------------
        % Section 6.2.5 Equation 36 ECMA-418-2:2025
        % ACF spectrum in the lag window [Phi'_z,tau(k)]
        magFFTlagWindowACF = abs(fft(lagWindowACF, 2*max(blockSize), 1));
        magFFTlagWindowACF(isnan(magFFTlagWindowACF)) = 0;

        % Section 6.2.5 Equation 37 ECMA-418-2:2025 [Nhat'_tonal(z)]
        % first estimation of specific partial loudness of tonal component in critical band
        bandTonalLoudness = meanScaledACF(1, :);
        mask = 2*max(magFFTlagWindowACF, [], 1)/(M/2) <= meanScaledACF(1, :);
        bandTonalLoudness(mask) = 2*max(magFFTlagWindowACF(:, mask), [], 1)/(M/2);

        % Section 6.2.5 Equation 38 & 39 ECMA-418-2:2025
        % [k_max(z)]
        [~, kz_max] = max(magFFTlagWindowACF, [], 1);
        % frequency of maximum tonal component in critical band [f_ton(z)]
        bandTonalFreqs = (kz_max - 1)*(sampleRate48k/(2*max(blockSize)));

        % Section 6.2.7 Equation 41 ECMA-418-2:2025 [N'_signal(l,z)]
        % specific partial loudness of complete band-pass signal in critical band
        bandLoudness = meanScaledACF(1, :);

        % Resampling to common time basis Section 6.2.6 ECMA-418-2:2025
        if i_interp(zBand) > 1
            % Note: use of interpolation function avoids rippling caused by
            % resample function, which otherwise affects specific loudness 
            % calculation for tonal and noise components
            l_n = size(meanScaledACF, 2);
            x = linspace(1, l_n, l_n);
            xq = linspace(1, l_n, i_interp(zBand)*(l_n - 1) + 1);

            bandTonalLoudness = interp1(x, bandTonalLoudness, xq);
            bandLoudness = interp1(x, bandLoudness, xq);
            bandTonalFreqs = interp1(x, bandTonalFreqs, xq);

        end  % end of if branch for interpolation

        % Remove end zero-padded samples Section 6.2.6 ECMA-418-2:2025
        l_end = ceil(size(pTarget_re, 1)/sampleRate48k*sampleRate1875) + 1;  % Equation 40 ECMA-418-2:2025
        bandTonalLoudness = bandTonalLoudness(1:l_end);
        bandLoudness = bandLoudness(1:l_end);
        bandTonalFreqs = bandTonalFreqs(1:l_end);

        % Noise reduction Section 6.2.7 ECMA-418-2:2020
        % ---------------------------------------------
        % Equation 42 ECMA-418-2:2025 signal-noise-ratio first approximation
        % (ratio of tonal component partial loudness to non-tonal component partial loudness in critical band)
        % [SNRhat(l,z)]
        SNRlz1 = bandTonalLoudness./((bandLoudness - bandTonalLoudness) + epsilon);

        % Equation 43 ECMA-418-2:2025 low pass filtered specific loudness
        % of non-tonal component in critical band [Ntilde'_tonal(l,z)]
        bandTonalLoudness = shmNoiseRedLowPass(bandTonalLoudness, sampleRate1875);

        % Equation 44 ECMA-418-2:2025 lowpass filtered SNR (improved estimation)
        % [SNRtilde(l,z)]
        SNRlz = shmNoiseRedLowPass(SNRlz1, sampleRate1875);

        % Equation 46 ECMA-418-2:2025 [g(z)]
        gz = csz_b(zBand)/(bandCentreFreqs(zBand)^dsz_b(zBand));

        % Equation 45 ECMA-418-2:2025 [nr(l,z)]
        crit = exp(-alpha*((SNRlz/gz) - beta));
        nrlz = 1 - crit;  % sigmoidal weighting function
        nrlz(crit >= 1) = 0;

        % Equation 47 ECMA-418-2:2025 [N'_tonal(l,z)]
        bandTonalLoudness = nrlz.*bandTonalLoudness;

        % Section 6.2.8 Equation 48 ECMA-418-2:2025 [N'_noise(l,z)]
        bandNoiseLoudness = shmNoiseRedLowPass(bandLoudness, sampleRate1875) - bandTonalLoudness;  % specific loudness of non-tonal component in critical band

        % Store critical band results
        % ---------------------------
        % specific time-dependent signal-noise-ratio in each critical band
%         specSNR(:, zBand, chan) = SNRlz1;

        % specific time-dependent partial loudness of signal in each critical band
%         specPartLoudness(:, zBand, chan) = bandLoudness;

        % specific time-dependent partial loudness of tonal component in each critical band  [N'_tonal(l,z)]
        specPartTonalLoudness(:, zBand, chan) = bandTonalLoudness;

        % specific time-dependent partial loudness of non-tonal component in each critical band [N'_noise(l,z)]
        specPartNoiseLoudness(:, zBand, chan) = bandNoiseLoudness;

        % time-dependent frequency of tonal component in each critical band [f_ton(z)]
        specPartTonalityFreqs(:, zBand, chan) = bandTonalFreqs;

    end  % end of for loop over ACF bands

    % set any tiny negative loudness values to 0
    specPartTonalLoudness(specPartTonalLoudness < 0) = 0;
    specPartNoiseLoudness(specPartNoiseLoudness < 0) = 0;

    % Calculation of specific partial tonality
    % ----------------------------------------
    % Section 6.2.8 Equation 49 ECMA-418-2:2025 [SNR(l)]
    overallSNR = max(specPartTonalLoudness, [], 2)./(sum(specPartNoiseLoudness, 2) + epsilon);  % loudness signal-noise-ratio

    % Section 6.2.8 Equation 50 ECMA-418-2:2025 [q(l)]
    crit = exp(-A*(overallSNR - B));
    ql = 1 - crit;  % sigmoidal scaling factor
    ql(crit >= 1) = 0;

    % Section 6.2.8 Equation 51 ECMA-418-2:2025 [T'(l,z)]
    specPartTonality = cal_T*cal_Tx*ql.*specPartTonalLoudness;  % time-dependent specific partial tonality

    % Calculation of time-averaged specific partial tonality Section 6.2.9
    % ECMA-418-2:2025 [T'(z)]
    for zBand = 53:-1:1
        mask = specPartTonality(:, zBand, chan) > 0.02;  % criterion Section 6.2.9 point 2
        mask(1:(58 - 1)) = 0;  % criterion Section 6.2.9 point 1

        % Section 6.2.9 Equation 53 ECMA-418-2:2025
        specPartTonalityAvg(1, zBand, chan)...
            = sum(specPartTonality(mask, zBand, chan), 1)./(nnz(mask) + epsilon);
        specPartTonalityAvgFreqs(1, zBand, chan)...
            = sum(specPartTonalityFreqs(mask, zBand, chan), 1)./(nnz(mask) + epsilon);
    end

    % Calculation of overall partial tonality Section 6.2.10
    % ----------------------------------------------
    % Further update can add the user input frequency range to determine
    % total tonality - not yet incorporated

    % Section 6.2.8 Equation 52 ECMA-418-2:2025
    % time (s) corresponding with results output [t]
    timeOut = transpose((0:(size(specPartTonality, 1) - 1))/sampleRate1875);

    % Section 6.2.10 Equation 61 ECMA-418-2:2025
    % Time-dependent total partial tonality [T(l)]
    [partTonalityTDep(:, chan), zmax] = max(specPartTonality(:, :, chan),...
                                           [], 2);
    for ll = size(specPartTonalityFreqs, 1):-1:1
        partTonalityTDepFreqs(ll, chan) = specPartTonalityFreqs(ll, zmax(ll), chan);
    end

    % Calculation of representative values Section 6.2.11 ECMA-418-2:2025
    % Time-averaged total partial tonality
    mask = partTonalityTDep(:, chan) > 0.02;  % criterion Section 6.2.9 point 2
    mask(1:(58 - 1)) = 0;    % criterion Section 6.2.9 point 1

    % Section 6.2.11 Equation 63 ECMA-418-2:2025
    % Time-averaged total partial tonality [T]
    partTonalityAvg(chan) = sum(partTonalityTDep(mask, chan))/(nnz(mask) + epsilon);

    if annoyWeight
        % Tonal annoyance weighted results
        annoyWeightx = ones(size(bandCentreFreqs));
        annoyWeightx(bandCentreFreqs...
                     >= 1e3) = 2.3*(log10(bandCentreFreqs(bandCentreFreqs...
                                                          >= 1e3)) - 3) + 1;
        specPartTonalityAnnoy(:, :, chan) = specPartTonality(:, :, chan).*annoyWeightx;

        for zBand = 53:-1:1
            mask = specPartTonalityAnnoy(:, zBand, chan) > 0.02;
            mask(1:(58 - 1)) = 0;
            specPartTonalityAvgAnnoy(1, zBand, chan)...
                = sum(specPartTonalityAnnoy(mask, zBand, chan), 1)./(nnz(mask) + epsilon);
        end

        partTonalityTDepAnnoy(:, chan) = max(specPartTonalityAnnoy(:, :, chan), [], 2);
        mask = partTonalityTDepAnnoy(:, chan) > 0.02;
        mask(1:(58 - 1)) = 0;
        partTonalityAvgAnnoy(chan) = sum(partTonalityTDepAnnoy(mask, chan))/(nnz(mask) + epsilon);
    end

    if waitBar
        close(w)  % close waitbar
    end

%% Output plotting

    % Plot figures
    % ------------
    if outPlot
        cmap_plasma = load('cmap_plasma.txt');
        % Plot results
        chan_lab = chans(chan);
        fig = figure;
        tiledlayout(fig, 2, 1);
        movegui(fig, 'center');
        ax1 = nexttile(1);
        surf(ax1, timeOut, bandCentreFreqs, permute(specPartTonality(:, :, chan),...
                                              [2, 1, 3]),...
             'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        ax1.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        ax1.YLim = [bandCentreFreqs(1), bandCentreFreqs(end)];
        ax1.CLim = [0, ceil(max(partTonalityTDep(:, chan))*10)/10];
        ax1.YTick = [63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3]; 
        ax1.YTickLabel = ["63", "125", "250", "500", "1k", "2k", "4k",...
                          "8k", "16k"]; 
        ax1.YScale = 'log';
        ax1.YLabel.String = "Frequency, Hz";
        ax1.XLabel.String = "Time, s";
        ax1.FontName =  'Arial';
        ax1.FontSize = 12;
        colormap(cmap_plasma);
        h = colorbar;
        set(get(h,'label'),'string', {'Specific partial tonality,'; 'tu_{SHM}/Bark_{SHM}'});

        % Create A-weighting filter
        weightFilt = weightingFilter('A-weighting', sampleRate48k);
        % Filter signal to determine A-weighted time-averaged level
        pA = weightFilt(pTarget_re(:, chan));
        LAeq = 20*log10(rms(pA)/2e-5);
        title(strcat(chan_lab,...
                     " signal sound pressure level =", {' '},...
                     num2str(round(LAeq,1)), "dB {\itL}_{Aeq}"),...
                     'FontWeight', 'normal', 'FontName', 'Arial');

        ax2 = nexttile(2);
        plot(ax2, timeOut, partTonalityAvg(1, chan)*ones(size(timeOut)), 'color', cmap_plasma(34, :),...
             'LineWidth', 1, 'DisplayName', "Time-" + string(newline) + "average");
        hold on
        plot(ax2, timeOut, partTonalityTDep(:, chan), 'color',  cmap_plasma(166, :),...
             'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold off
        ax2.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        ax2.YLim = [0, 1.1*ceil(max(partTonalityTDep(:, chan))*10)/10];
        ax2.XLabel.String = "Time, s";
        ax2.YLabel.String = "Partial tonality, tu_{SHM}";
        ax2.XGrid = 'on';
        ax2.YGrid = 'on';
        ax2.GridAlpha = 0.075;
        ax2.GridLineStyle = '--';
        ax2.GridLineWidth = 0.25;
        ax2.FontName = 'Arial';
        ax2.FontSize = 12;
        lgd = legend('Location', 'eastoutside', 'FontSize', 8);
        lgd.Title.String = "Overall";

        if annoyWeight
            % Plot tonal annoyance-weighted results
            chan_lab = chans(chan);
            fig = figure;
            tiledlayout(fig, 2, 1);
            movegui(fig, 'center');
            ax1 = nexttile(1);
            surf(ax1, timeOut, bandCentreFreqs, permute(specPartTonalityAnnoy(:, :, chan),...
                                                  [2, 1, 3]),...
                 'EdgeColor', 'none', 'FaceColor', 'interp');
            view(2);
            ax1.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
            ax1.YLim = [bandCentreFreqs(1), bandCentreFreqs(end)];
            ax1.CLim = [0, ceil(max(partTonalityTDep(:, chan))*10)/10];
            ax1.YTick = [63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3]; 
            ax1.YTickLabel = ["63", "125", "250", "500", "1k", "2k", "4k",...
                              "8k", "16k"]; 
            ax1.YScale = 'log';
            ax1.YLabel.String = "Frequency, Hz";
            ax1.XLabel.String = "Time, s";
            ax1.FontName =  'Arial';
            ax1.FontSize = 12;
            colormap(cmap_plasma);
            h = colorbar;
            set(get(h,'label'),'string', {'Specific partial tonality'; '(tonal annoyance-weighted),'; 'tu_{SHMaw}/Bark_{SHM}'});
    
            % Create A-weighting filter
            weightFilt = weightingFilter('A-weighting', sampleRate48k);
            % Filter signal to determine A-weighted time-averaged level
            pA = weightFilt(pTarget_re(:, chan));
            LAeq = 20*log10(rms(pA)/2e-5);
            title(strcat(chan_lab,...
                         " signal sound pressure level =", {' '},...
                         num2str(round(LAeq,1)), "dB {\itL}_{Aeq}"),...
                         'FontWeight', 'normal', 'FontName', 'Arial');
    
            ax2 = nexttile(2);
            plot(ax2, timeOut, partTonalityAvgAnnoy(1, chan)*ones(size(timeOut)),...
                 ':', 'color', cmap_plasma(34, :),...
                 'LineWidth', 1.5, 'DisplayName', "Time-" + string(newline) + "average");
            hold on
            plot(ax2, timeOut, partTonalityTDepAnnoy(:, chan), 'color',  cmap_plasma(166, :),...
                 'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
            hold off
            ax2.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
            ax2.YLim = [0, 1.1*ceil(max(partTonalityTDep(:, chan))*10)/10];
            ax2.XLabel.String = "Time, s";
            ax2.YLabel.String = "Partial tonality (tonal annoyance-weighted), tu_{SHMaw}";
            ax2.XGrid = 'on';
            ax2.YGrid = 'on';
            ax2.GridAlpha = 0.075;
            ax2.GridLineStyle = '--';
            ax2.GridLineWidth = 0.25;
            ax2.FontName = 'Arial';
            ax2.FontSize = 12;
            lgd = legend('Location', 'eastoutside', 'FontSize', 8);
            lgd.Title.String = "Overall";
        end
    end

end  % end of for loop over channels

% Discard singleton dimensions (must come after channels loop)
if chansIn > 1
    specPartTonalityAvg = squeeze(specPartTonalityAvg);
    specPartTonalityAvgFreqs = squeeze(specPartTonalityAvgFreqs);
    if annoyWeight
        specPartTonalityAvgAnnoy = squeeze(specPartTonalityAvgAnnoy);
    end
else
    specPartTonalityAvg = transpose(specPartTonalityAvg);
    specPartTonalityAvgFreqs = transpose(specPartTonalityAvgFreqs);
    if annoyWeight
        specPartTonalityAvgAnnoy = transpose(specPartTonalityAvgAnnoy);
    end
end

%% Output assignment

% Assign outputs to structure
partTonalitySHM.specPartTonality = specPartTonality;
partTonalitySHM.specPartTonalityAvg = specPartTonalityAvg;
partTonalitySHM.specPartTonalityFreqs = specPartTonalityFreqs;
partTonalitySHM.specPartTonalityAvgFreqs = specPartTonalityAvgFreqs;
partTonalitySHM.specPartTonalLoudness = specPartTonalLoudness;
partTonalitySHM.specNoiseLoudness = specPartNoiseLoudness;
partTonalitySHM.partTonalityTDep = partTonalityTDep;
partTonalitySHM.partTonalityAvg = partTonalityAvg;
partTonalitySHM.partTonalityTDepFreqs = partTonalityTDepFreqs;
partTonalitySHM.bandCentreFreqs = bandCentreFreqs;
partTonalitySHM.timeOut = timeOut;
partTonalitySHM.soundField = soundField;

% additional outputs for tonal annoyance weighting
if annoyWeight
    partTonalitySHM.specPartTonalityAnnoy = specPartTonalityAnnoy;
    partTonalitySHM.specPartTonalityAvgAnnoy = specPartTonalityAvgAnnoy;
    partTonalitySHM.partTonalityTDepAnnoy = partTonalityTDepAnnoy;
    partTonalitySHM.partTonalityAvgAnnoy = partTonalityAvgAnnoy;
end

% end of function
