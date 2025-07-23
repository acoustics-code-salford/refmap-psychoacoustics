function tonalitySHM = acousticSHMTonality(p, sampleRateIn, axisN, soundField, waitBar, outPlot)
% tonalitySHM = acousticSHMTonality(p, sampleRateIn, axisN, soundField, waitBar, outPlot)
%
% Returns tonality values and frequencies according to ECMA-418-2:2025
% (using the Sottek Hearing Model) for an input calibrated single mono
% or single stereo audio (sound pressure) time-series signal, p.
%
% Inputs
% ------
% p : vector or 2D matrix
%     the input signal as single mono or stereo audio (sound
%     pressure) signals
%
% sampleRateIn : integer
%                the sample rate (frequency) of the input signal(s)
%
% axisN : integer (1 or 2, default: 1)
%         the time axis along which to calculate the tonality
%
% soundField : keyword string (default: 'freeFrontal')
%              determines whether the 'freeFrontal' or 'diffuse' field stages
%              are applied in the outer-middle ear filter, or 'noOuter' uses
%              only the middle ear stage, or 'noEar' omits ear filtering.
%              Note: these last two options are beyond the scope of the
%              standard, but may be useful if recordings made using
%              artificial outer/middle ear are to be processed using the
%              specific recorded responses.
%
% waitBar : keyword string (default: true)
%           determines whether a progress bar displays during processing
%           (set waitBar to false for doing multi-file parallel calculations)
%
% outPlot : Boolean true/false (default: false)
%           flag indicating whether to generate a figure from the output
%           (set outPlot to false for doing multi-file parallel calculations)
% 
% Returns
% -------
%
% tonalitySHM : structure
%               contains the output
%
% tonalitySHM contains the following outputs:
%
% specTonality : matrix
%                time-dependent specific tonality for each (half) critical
%                band
%                arranged as [time, bands(, channels)]
%
% specTonalityFreqs : matrix
%                     time-dependent frequencies of the dominant tonal
%                     components corresponding with each of the
%                     time-dependent specific tonality values in each
%                     (half) critical band
%                     arranged as [time, bands(, channels)]
%
% specTonalityAvg : matrix
%                   time-averaged specific tonality for each (half)
%                   critical band
%                   arranged as [bands(, channels)]
%
% specTonalityAvgFreqs : matrix
%                        frequencies of the dominant tonal components
%                        corresponding with each of the
%                        time-averaged specific tonality values in each
%                        (half) critical band
%                        arranged as [bands(, channels)]
%
% specTonalLoudness : matrix
%                     time-dependent specific tonal loudness for each
%                     (half) critical band
%                     arranged as [time, bands(, channels)]
%
% specNoiseLoudness : matrix
%                     time-dependent specific noise loudness for each
%                     (half) critical band
%                     arranged as [time, bands(, channels)]
%
% tonalityTDep : vector or matrix
%                time-dependent overall tonality
%                arranged as [time(, channels)]
%
% tonalityTDepFreqs : vector or matrix
%                     time-dependent frequencies of the dominant tonal
%                     components corresponding with the
%                     time-dependent overall tonality values
%                     arranged as [time(, channels)]
%
% tonalityAvg : number or vector
%               time-averaged overall tonality
%               arranged as [tonality(, channels)]
%
% bandCentreFreqs : vector
%                   centre frequencies corresponding with each (half)
%                   critical band rate scale width
%
% timeOut : vector
%           time (seconds) corresponding with time-dependent outputs
%
% soundField : string
%              identifies the soundfield type applied (the input argument
%              soundField)
%
%
% If outPlot=true, a set of plots is returned illustrating the energy
% time-averaged A-weighted sound level, the time-dependent specific and
% overall tonality, with the latter also indicating the time-aggregated
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
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk) &
%          Matt Torjussen (matt@anv.co.uk)
% Institution: University of Salford / ANV Measurement Systems
%
% Date created: 07/08/2023
% Date last modified: 23/07/2025
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
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Input checks
% Orient input matrix
if axisN == 2
    p = p.';
end

% Check the length of the input data (must be longer than 300 ms)
if size(p, 1) <=  300/1000*sampleRateIn
    error('Error: Input signal is too short along the specified axis to calculate tonality (must be longer than 300 ms)')
end

% Check the channel number of the input data
if size(p, 2) > 2
    error('Error: Input signal comprises more than two channels')
else
    chansIn = size(p, 2);
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
    [p_re, ~] = shmResample(p, sampleRateIn);
else  % don't resample
    p_re = p;
end

% Section 5.1.2 ECMA-418-2:2025 Fade in weighting and zero-padding
pn = shmPreProc(p_re, max(blockSize), max(hopSize));

% Apply outer & middle ear filter
% -------------------------------
%
% Section 5.1.3.2 ECMA-418-2:2025 Outer and middle/inner ear signal filtering
pn_om = shmOutMidEarFilter(pn, soundField);

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
    pn_omz = shmAuditoryFiltBank(pn_om(:, chan));

    % Autocorrelation function analysis
    % ---------------------------------
    % Duplicate Banded Data for ACF
    % Averaging occurs over neighbouring bands, to do this the segmentation
    % needs to be duplicated for neigbouring bands. 'Dupe' has been added
    % to variables to indicate that the vectors/matrices have been modified
    % for duplicated neigbouring bands.
    
    pn_omzDupe = [pn_omz(:, 1:5), pn_omz(:, 2:18), pn_omz(:, 16:26),...
                   pn_omz(:, 26:53)];
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
        [pn_lz, ~] = shmSignalSegment(pn_omzDupe(:, zBand), 1,...
                                      blockSizeDupe(zBand), overlap, i_start);
 
        % Transformation into Loudness
        % ----------------------------
        % Sections 5.1.6 to 5.1.9 ECMA-418-2:2025 [N'_basis(z)]
        [pn_rlz, bandBasisLoudness, ~]...
            = shmBasisLoudness(pn_lz, bandCentreFreqsDupe(zBand));
        basisLoudnessArray{zBand} = bandBasisLoudness;

        % Apply ACF
        % ACF implementation using DFT
        % Section 6.2.2 Equations 27 & 28 ECMA-418-2:2025
        % [phi_unscaled,l,z(m)]
        unscaledACF = ifft(abs(fft(pn_rlz, 2*blockSizeDupe(zBand), 1)).^2,...
                           2*blockSizeDupe(zBand), 1);
        % Section 6.2.2 Equation 29 ECMA-418-2:2025 [phi_l,z(m)]
        denom = sqrt(cumsum(pn_rlz.^2, 1, 'reverse').*flipud(cumsum(pn_rlz.^2, 1)))...
                + epsilon;

        % note that the block length is used here, rather than the 2*s_b,
        % for compatability with the remaining code - beyond 0.75*s_b is
        % assigned (unused) zeros in the next line
        unbiasedNormACF = unscaledACF(1:blockSizeDupe(zBand), :)./denom;
        unbiasedNormACF((0.75*blockSizeDupe(zBand) + 1):blockSizeDupe(zBand), :) = 0;

        % Section 6.2.2 Equation 30 ECMA-418-2:2025 [phi_z'(m)
        unbiasedNormACFDupe{zBand} = basisLoudnessArray{zBand}.*unbiasedNormACF;

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

        % Estimation of tonal loudness
        % ----------------------------
        % Section 6.2.5 Equation 36 ECMA-418-2:2025
        % ACF spectrum in the lag window [Phi'_z,tau(k)]
        magFFTlagWindowACF = abs(fft(lagWindowACF, 2*max(blockSize), 1));
        magFFTlagWindowACF(isnan(magFFTlagWindowACF)) = 0;

        % Section 6.2.5 Equation 37 ECMA-418-2:2025 [Nhat'_tonal(z)]
        % first estimation of specific loudness of tonal component in critical band
        bandTonalLoudness = meanScaledACF(1, :);
        mask = 2*max(magFFTlagWindowACF, [], 1)/(M/2) <= meanScaledACF(1, :);
        bandTonalLoudness(mask) = 2*max(magFFTlagWindowACF(:, mask), [], 1)/(M/2);

        % Section 6.2.5 Equation 38 & 39 ECMA-418-2:2025
        % [k_max(z)]
        [~, kz_max] = max(magFFTlagWindowACF, [], 1);
        % frequency of maximum tonal component in critical band [f_ton(z)]
        bandTonalFreqs = (kz_max - 1)*(sampleRate48k/(2*max(blockSize)));

        % Section 6.2.7 Equation 41 ECMA-418-2:2025 [N'_signal(l,z)]
        % specific loudness of complete band-pass signal in critical band
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
        l_end = ceil(size(p_re, 1)/sampleRate48k*sampleRate1875) + 1;  % Equation 40 ECMA-418-2:2025
        bandTonalLoudness = bandTonalLoudness(1:l_end);
        bandLoudness = bandLoudness(1:l_end);
        bandTonalFreqs = bandTonalFreqs(1:l_end);

        % Noise reduction Section 6.2.7 ECMA-418-2:2020
        % ---------------------------------------------
        % Equation 42 ECMA-418-2:2025 signal-noise-ratio first approximation
        % (ratio of tonal component loudness to non-tonal component loudness in critical band)
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

        % specific time-dependent loudness of signal in each critical band
%         specLoudness(:, zBand, chan) = bandLoudness;

        % specific time-dependent loudness of tonal component in each critical band  [N'_tonal(l,z)]
        specTonalLoudness(:, zBand, chan) = bandTonalLoudness;

        % specific time-dependent loudness of non-tonal component in each critical band [N'_noise(l,z)]
        specNoiseLoudness(:, zBand, chan) = bandNoiseLoudness;

        % time-dependent frequency of tonal component in each critical band [f_ton(z)]
        specTonalityFreqs(:, zBand, chan) = bandTonalFreqs;

    end

    % Calculation of specific tonality
    % --------------------------------
    % Section 6.2.8 Equation 49 ECMA-418-2:2025 [SNR(l)]
    overallSNR = max(specTonalLoudness, [], 2)./(sum(specNoiseLoudness, 2) + epsilon);  % loudness signal-noise-ratio
    
    % Section 6.2.8 Equation 50 ECMA-418-2:2025 [q(l)]
    crit = exp(-A*(overallSNR - B));
    ql = 1 - crit;  % sigmoidal scaling factor
    ql(crit >= 1) = 0;
    
    % Section 6.2.8 Equation 51 ECMA-418-2:2025 [T'(l,z)]
    specTonality = cal_T*cal_Tx*ql.*specTonalLoudness;  % time-dependent specific tonality
    
    % Calculation of time-averaged specific tonality Section 6.2.9
    % ECMA-418-2:2025 [T'(z)]
    for zBand = 53:-1:1
        mask = specTonality(:, zBand, chan) > 0.02;  % criterion Section 6.2.9 point 2
        mask(1:(58 - 1)) = 0;  % criterion Section 6.2.9 point 1

        % Section 6.2.9 Equation 53 ECMA-418-2:2025
        specTonalityAvg(1, zBand, chan)...
            = sum(specTonality(mask, zBand, chan), 1)./(nnz(mask) + epsilon);
        specTonalityAvgFreqs(1, zBand, chan)...
            = sum(specTonalityFreqs(mask, zBand, chan), 1)./(nnz(mask) + epsilon);
    end

    % Calculation of overall tonality Section 6.2.10
    % ----------------------------------------------
    % Further update can add the user input frequency range to determine
    % total tonality - not yet incorporated

    % Section 6.2.8 Equation 52 ECMA-418-2:2025
    % time (s) corresponding with results output [t]
    timeOut = (0:(size(specTonality, 1) - 1))/sampleRate1875;

    % Section 6.2.10 Equation 61 ECMA-418-2:2025
    % Time-dependent total tonality [T(l)]
    [tonalityTDep(:, chan), zmax] = max(specTonality(:, :, chan),...
                                           [], 2);
    for ll = size(specTonalityFreqs, 1):-1:1
        tonalityTDepFreqs(ll, chan) = specTonalityFreqs(ll, zmax(ll), chan);
    end
    
    % Calculation of representative values Section 6.2.11 ECMA-418-2:2025
    % Time-averaged total tonality
    mask = tonalityTDep(:, chan) > 0.02;  % criterion Section 6.2.9 point 2
    mask(1:(58 - 1)) = 0;    % criterion Section 6.2.9 point 1

    % Section 6.2.11 Equation 63 ECMA-418-2:2025
    % Time-averaged total tonality [T]
    tonalityAvg(chan) = sum(tonalityTDep(mask, chan))/(nnz(mask) + epsilon);

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
        surf(ax1, timeOut, bandCentreFreqs, permute(specTonality(:, :, chan),...
                                              [2, 1, 3]),...
             'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        ax1.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        ax1.YLim = [bandCentreFreqs(1), bandCentreFreqs(end)];
        ax1.CLim = [0, ceil(max(tonalityTDep(:, chan))*10)/10];
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
        set(get(h,'label'),'string', {'Specific Tonality,'; 'tu_{SHM}/Bark_{SHM}'});

        % Create A-weighting filter
        weightFilt = weightingFilter('A-weighting', sampleRate48k);
        % Filter signal to determine A-weighted time-averaged level
        pA = weightFilt(p_re(:, chan));
        LAeq = 20*log10(rms(pA)/2e-5);
        title(strcat(chan_lab,...
                     " signal sound pressure level =", {' '},...
                     num2str(round(LAeq,1)), "dB {\itL}_{Aeq}"),...
                     'FontWeight', 'normal', 'FontName', 'Arial');
        
        ax2 = nexttile(2);
        plot(ax2, timeOut, tonalityAvg(1, chan)*ones(size(timeOut)), 'color', cmap_plasma(34, :),...
             'LineWidth', 1, 'DisplayName', "Time-" + string(newline) + "average");
        hold on
        plot(ax2, timeOut, tonalityTDep(:, chan), 'color',  cmap_plasma(166, :),...
             'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold off
        ax2.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        ax2.YLim = [0, 1.1*ceil(max(tonalityTDep(:, chan))*10)/10];
        ax2.XLabel.String = "Time, s";
        ax2.YLabel.String = "Tonality, tu_{SHM}";
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

end  % end of for loop over channels

%% Output assignment

% Discard singleton dimensions
if chansIn > 1
    specTonalityAvg = squeeze(specTonalityAvg);
    specTonalityAvgFreqs = squeeze(specTonalityAvgFreqs);
else
    specTonalityAvg = transpose(specTonalityAvg);
    specTonalityAvgFreqs = transpose(specTonalityAvgFreqs);
end

% Assign outputs to structure
tonalitySHM.specTonality = specTonality;
tonalitySHM.specTonalityAvg = specTonalityAvg;
tonalitySHM.specTonalityFreqs = specTonalityFreqs;
tonalitySHM.specTonalityAvgFreqs = specTonalityAvgFreqs;
tonalitySHM.specTonalLoudness = specTonalLoudness;
tonalitySHM.specNoiseLoudness = specNoiseLoudness;
tonalitySHM.tonalityTDep = tonalityTDep;
tonalitySHM.tonalityAvg = tonalityAvg;
tonalitySHM.tonalityTDepFreqs = tonalityTDepFreqs;
tonalitySHM.bandCentreFreqs = bandCentreFreqs;
tonalitySHM.timeOut = timeOut;
tonalitySHM.soundField = soundField;

% end of function
