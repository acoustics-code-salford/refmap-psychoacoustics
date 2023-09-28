function [tonalityAvg, tonalityTimeVar, tonalityTimeVarFreqs,...
          specificTonality, specificTonalityFreqs, specificTonalityAvg,...
          specificTonalityAvgFreqs, bandCentreFreqs,...
          specificTonalLoudness, specificNoiseLoudness]...
          = acousticHMSTonality_(p, sampleRatein, axisn, outplot, ecma)
% [tonalityAvg, tonalityTimeVar, tonalityTimeVarFreqs, specificTonality,
%  specificTonalityFreqs, specificTonalityAvg, specificTonalityAvgFreqs,
%  bandCentreFreqs]
%  = acousticHMSTonality(p, sampleRatein, axisn, outplot, ecma)
%
% Returns tonality values and frequencies according to ECMA-418-2:2022
% (using the Hearing Model of Sottek) for an input calibrated single mono
% or single stereo audio (sound pressure) time-series signal, p.
%
% Inputs
% ------
% p : vector or 2D matrix
%     the input signal as single mono or stereo audio (sound
%     pressure) signals
%
% sampleRatein : integer
%                the sample rate (frequency) of the input signal(s)
%
% axisn : integer (1 or 2, default: 1)
%         the time axis along which to calculate the tonality
%
% outplot : Boolean true/false (default: false)
%           flag indicating whether to generate a figure from the output
%
% ecma : Boolean true/false (default: true)
%        flag indicating whether to maintain strict standard adherence to
%        ECMA-418-2:2022 Equation 40, or otherwise to use an alternative
%        that provides closer time-alignment of the time-dependent tonality
%        with the original signal
% 
% Returns
% -------
% For each channel in the input signal:
%
% tonalityAvg : number or vector
%               average (overall) tonality value
% 
% tonalityTimeVar : vector or 2D matrix
%                   time-dependent overall tonality values
%
% tonalityTimeVarFreqs : vector or 2D matrix
%                        time-dependent frequencies of the dominant tonal
%                        components corresponding with the time-dependent
%                        overall tonality values
%
% specificTonality : 2D or 3D matrix
%                    time-dependent specific tonality values in each
%                    half-critical band rate scale width
%
% specificTonalityFreqs : 2D or 3D matrix
%                         time-dependent frequencies of the dominant tonal
%                         components corresponding with each of the
%                         time-dependent specific tonality values in each
%                         half-critical band rate scale width
%
% specificTonalityAvg : vector or 2D matrix
%                       time-averaged specific tonality values in each
%                       half-critical band rate scale width
%
% specificTonalityAvgFreqs : vector or 2D matrix
%                            time-averaged frequencies of the dominant
%                            tonal components corresponding with each of
%                            the half-critical band rate scale width
%
% bandCentreFreqs : vector
%                   centre frequencies corresponding with each half-Bark
%                   critical band rate scale width
%
% specificTonalLoudness : 2D or 3D matrix
%                         time-dependent specific loudness of the tonal
%                         components in each half-critical band rate scale
%                         width
%
% specificNoiseLoudness : 2D or 3D matrix
%                         time-dependent specific loudness of the noise
%                         components in each half-critical band rate scale
%                         width
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
% Date last modified: 27/09/2023
% MATLAB version: 2022b
%
% Copyright statement: This file and code is part of work undertaken within
% the RefMap project (www.refmap.eu), and is subject to licence as detailed
% in the code repository
% (https://github.com/acoustics-code-salford/refmap-psychoacoustics)
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
        sampleRatein (1, 1) double {mustBePositive, mustBeInteger}
        axisn (1, 1) {mustBeInteger, mustBeInRange(axisn, 1, 2)} = 1
        outplot {mustBeNumericOrLogical} = false
        ecma {mustBeNumericOrLogical} = true
    end

%% Load path
addpath(genpath("refmap-psychoacoustics\src\mlab\"))

%% Input checks
% Orient input matrix
if axisn == 2
    p = p.';
end

% Check the length of the input data (must be longer than 300 ms)
if size(p, 1) <=  300/1000*sampleRatein
    error('Error: Input signal is too short to calculate tonality (must be longer than 300 ms)')
end

% Check the channel number of the input data
if size(p, 2) > 2
    error('Error: Input signal comprises more than two channels')
else
    inchans = size(p, 2);
    if inchans > 1
        chans = ["Stereo left";
                 "Stereo right"];
    else
        chans = "Mono";
    end
end

%% Define constants

sampleRate48k = 48e3;  % Signal sample rate prescribed to be 48kHz (to be used for resampling), Section 5.1.1 ECMA-418-2:2022
deltaFreq0 = 81.9289;  % defined in Section 5.1.4.1 ECMA-418-2:2022
c = 0.1618;  % Half-Bark band centre-frequency demoninator constant defined in Section 5.1.4.1 ECMA-418-2:2022

halfBark = 0.5:0.5:26.5;  % half-critical band rate scale
bandCentreFreqs = (deltaFreq0/c)*sinh(c*halfBark);  % Section 5.1.4.1 Equation 9 ECMA-418-2:2022
dfz = sqrt(deltaFreq0^2 + (c*bandCentreFreqs).^2);  % Section 5.1.4.1 Equation 10 ECMA-418-2:2022

% Block and hop sizes Section 6.2.2 Table 4 ECMA-418-2:2022
overlap = 0.75;  % block overlap proportion
% block sizes
blockSize = [8192*ones(1, 3), 4096*ones(1, 13), 2048*ones(1, 9), 1024*ones(1, 28)];
hopSize = (1 - overlap)*blockSize;  % hop sizes (section 5.1.2 footnote 3 ECMA 418-2:2022)

sampleRate1875 = sampleRate48k/min(hopSize);  % Output sample rate based on hop sizes - Resampling to common time basis Section 6.2.6 ECMA-418-2:2022

% Number of bands that need averaging. Section 6.2.3 Table 5 ECMA-418-2:2022
NBandsAvg = [0, 1, 2*ones(1,14), ones(1,9), zeros(1,28);...
             1, 1, 2*ones(1,14), ones(1,9), zeros(1,28)];

% Critical band interpolation factors from Section 6.2.6 Table 6 ECMA-418-2:2022
i_interp = blockSize/min(blockSize);

% Noise reduction constants from Section 6.2.7 Table 7 ECMA-418-2:2022
alpha = 20;
beta = 0.07;

% Sigmoid function factor parameters Section 6.2.7 Table 8 ECMA-418-2:2022
csz_b = [18.21*ones(1, 3), 12.14*ones(1, 13), 417.54*ones(1, 9),...
         962.68*ones(1, 28)]; 
dsz_b = [0.36*ones(1, 3), 0.36*ones(1, 13), 0.71*ones(1, 9),...
         0.69*ones(1, 28)]; 

% Scaling factor constants from Section 6.2.8 Table 9 ECMA-418-2:2022
A = 35;
B = 0.003;
cal_T = 2.8785151;  % calibration factor in Section 6.2.8 Equation 51 ECMA-418-2:2022
cal_Tx = 1/1.00186;  % Adjustment to calibration factor cal_T Footnote 22 ECMA-418-2:2022

%% Signal processing

% Input pre-processing
% --------------------
if sampleRatein ~= sampleRate48k  % Resample signal
    [p_re, ~] = acousticHMSResample_(p, sampleRatein);
else  % don't resample
    p_re = p;
end

% Section 5.1.2 ECMA-418-2:2022 Fade in weighting and zero-padding
pn = acousticHMSPreProc_(p_re, max(blockSize), max(hopSize));

% Apply outer & middle ear filter
% -------------------------------
%
% Section 5.1.3.2 ECMA-418-2:2022 Outer and middle/inner ear signal filtering
pn_om = acousticHMSOutMidEarFilter_(pn);

n_steps = 115;  % approximate number of calculation steps

% Loop through channels in file
% -----------------------------
for chan = size(pn_om, 2):-1:1
    w = waitbar(0, "Initialising...");
    i_step = 1;

    % Apply auditory filter bank
    % --------------------------
    waitbar(i_step/n_steps, w, 'Applying auditory filters...');
    i_step = i_step + 1;

    % Filter equalised signal using 53 1/2Bark ERB filters according to 
    % Section 5.1.4.2 ECMA-418-2:2022
    pn_omz = acousticHMSAuditoryFiltBank_(pn_om(:, chan), false);

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
    
        waitbar(((62 - zBand) + i_step)/n_steps, w, strcat("Applying ACF in 61 bands, ",...
            num2str(zBand), ' to go...'));
        
        % Segmentation into blocks
        % ------------------------
        % Section 5.1.5 ECMA-418-2:2022
        i_start = blockSizeDupe(1) - blockSizeDupe(zBand) + 1;
        pn_lz = signalSegment_(pn_omzDupe(:, zBand), 1,...
                               blockSizeDupe(zBand), overlap, i_start);
 
        % Transformation into Loudness
        % ----------------------------
        % Sections 5.1.6 to 5.1.9 ECMA-418-2:2022
        [pn_rlz, bandBasisLoudness, ~]...
            = acousticHMSBandBasisLoudness_(pn_lz, bandCentreFreqsDupe(zBand));
        basisLoudness{zBand} = bandBasisLoudness;

        % Apply ACF
        % ACF implementation using DFT
        % Section 6.2.2 Equations 27 & 28 ECMA-418-2:2022
        unscaledACF = ifft(abs(fft(pn_rlz, 2*blockSizeDupe(zBand), 1)).^2,...
                           2*blockSizeDupe(zBand), 1);
        % Section 6.2.2 Equation 29 ECMA-418-2:2022
        denom = sqrt(cumsum(pn_rlz.^2, 1, 'reverse').*flipud(cumsum(pn_rlz.^2)))...
                + 1e-12; 
        unbiasedNormACF = unscaledACF(1:blockSizeDupe(zBand), :)./denom;  % note that the block length is used here, rather than the 2*s_b, for compatability with the remaining code - beyond 0.75*s_b is assigned (unused) zeros in the next line
        unbiasedNormACF((0.75*blockSizeDupe(zBand) + 1):blockSizeDupe(zBand), :) = 0;

        % Section 6.2.2 Equation 30 ECMA-418-2:202
        unbiasedNormACFDupe{zBand} = basisLoudness{zBand}.*unbiasedNormACF;

    end
    
    i_step = i_step + 62;  % increment calculation step for waitbar

    % Average the ACF over nB bands - Section 6.2.3 ECMA-418-2:2022        
    for zBand = 53:-1:1  % Loop through 53 critical band filtered signals
        waitbar(((54 - zBand) + i_step)/n_steps, w,...
                strcat("Calculating sound quality in 53 bands, ",...
                       num2str(zBand), " to go..."));
        
        NBZ = NBandsAvg(1, zBand) + NBandsAvg(2, zBand) + 1; % Total number of bands to average over
        
        % Averaging of frequency bands
        meanScaledACF = mean(reshape(cell2mat(unbiasedNormACFDupe(i_NBandsAvgDupe(1, zBand):i_NBandsAvgDupe(2, zBand))),...
                                     blockSize(zBand), [], NBZ), 3);

        % Average the ACF over adjacent time blocks
        if zBand <= 16 
            meanScaledACF = movmean(meanScaledACF, 3, 2, 'omitnan',...
                                    'EndPoints', 'fill');
        end
        
        % Application of ACF lag window Section 6.2.4 ECMA-418-2:2022
        tauz_start = max(0.5/dfz(zBand), 2e-3);  % Equation 31 ECMA-418-2:2022
        tauz_end = max(4/dfz(zBand), tauz_start + 1e-3);  % Equation 32 ECMA-418-2:2022
        % Equations 33 & 34 ECMA-418-2:2022
        mz_start = ceil(tauz_start*sampleRate48k);  % Starting lag window index
        mz_end = floor(tauz_end*sampleRate48k);  % Ending lag window index
        M = mz_end - mz_start + 1;
        % Equation 35 ECMA-418-2:2022
        lagWindowACF = zeros(size(meanScaledACF));
        % lag-windowed, detrended ACF
        lagWindowACF(mz_start:mz_end, :) = meanScaledACF(mz_start:mz_end, :)...
                                           - mean(meanScaledACF(mz_start:mz_end, :));
        
        % Estimation of tonal loudness
        % ----------------------------
        % Section 6.2.5 Equation 36 ECMA-418-2:2022
        % ACF spectrum in the lag window
        magFFTlagWindowACF = abs(fft(lagWindowACF, 2*max(blockSize), 1));
        magFFTlagWindowACF(isnan(magFFTlagWindowACF)) = 0;
    
        % Section 6.2.5 Equation 37 ECMA-418-2:2022
        % first estimation of specific loudness of tonal component in critical band
        bandTonalLoudness = meanScaledACF(1, :);
        mask = 2*max(magFFTlagWindowACF, [], 1)/(M/2) <= meanScaledACF(1, :);
        bandTonalLoudness(mask) = 2*max(magFFTlagWindowACF(:, mask), [], 1)/(M/2);
    
        % Section 6.2.5 Equations 38 & 39 ECMA-418-2:2022
        [~, kz_max] = max(magFFTlagWindowACF, [], 1);
        % frequency of maximum tonal component in critical band
        bandTonalFreqs = kz_max*(sampleRate48k/(2*max(blockSize)));

        % Section 6.2.7 Equation 41 ECMA-418-2:2022
        bandLoudness = meanScaledACF(1, :);  % specific loudness of complete band-pass signal in critical band
        
        % Resampling to common time basis Section 6.2.6 ECMA-418-2:2022
        if i_interp(zBand) > 1
            % Note: use of interpolation function avoids rippling caused by
            % resample function, which disrupts specific loudness 
            % calculation for tonal and noise components
            l_n = size(meanScaledACF, 2);
            x = linspace(1, l_n, l_n);
            xq = linspace(1, l_n, i_interp(zBand)*l_n);
            bandTonalLoudness = interp1(x, bandTonalLoudness, xq);
            bandLoudness = interp1(x, bandLoudness, xq);
            bandTonalFreqs = interp1(x, bandTonalFreqs, xq);

        end

        % Remove end zero-padded samples Section 6.2.6 ECMA-418-2:2022
        % Note: in this part of the standard, the effect of the end
        % zero-padding to the input is removed from the processed signals.
        % There is no mention of realigning the processed signals to
        % compensate for the start zero-padding to the input. The
        % alternative terms below incorporate an amendment to the standard
        % to compensate for the start zero-padding, which results in
        % improved time alignment of the processed signals with the input.
        if ecma == true
            l_end = ceil(size(p, 1)/sampleRate48k*sampleRate1875) + 1;  % Equation 40 ECMA-418-2:2022
            l_start = 1;  % start block
        else
            l_end = ceil((size(p, 1) + max(blockSize))/sampleRate48k*sampleRate1875) + 1;  % Equation 40 ECMA-418-2:2022 (edited to account for start zero-padding)
            l_start = floor(max(blockSize)/sampleRate48k*sampleRate1875) + 1;  % Additional term to remove start zero-padding lag
        end

        bandTonalLoudness = bandTonalLoudness(l_start:l_end);
        bandLoudness = bandLoudness(l_start:l_end);
        bandTonalFreqs = bandTonalFreqs(l_start:l_end);

        % Noise reduction Section 6.2.7 ECMA-418-2:2020
        % ---------------------------------------------
        SNRlz1 = bandTonalLoudness./((bandLoudness - bandTonalLoudness) + 1e-12);  % Equation 42 ECMA-418-2:2022 signal-noise-ratio first approximation (ratio of tonal component loudness to non-tonal component loudness in critical band)
        bandTonalLoudness = acousticHMSNoiseRedLowPass_(bandTonalLoudness, sampleRate1875);  % Equation 43 ECMA-418-2:2022 low pass filtered specific loudness of non-tonal component in critical band
        SNRlz = acousticHMSNoiseRedLowPass_(SNRlz1, sampleRate1875);  % Equation 44 ECMA-418-2:2022 lowpass filtered SNR (improved estimation)
        gz = csz_b(zBand)/(bandCentreFreqs(zBand)^dsz_b(zBand));  % Equation 46 ECMA-418-2:2022
        % Equation 45 ECMA-418-2:2022
        crit = exp(-alpha*((SNRlz/gz) - beta));
        nrlz = 1 - crit;  % sigmoidal weighting function
        nrlz(crit >= 1) = 0;
        bandTonalLoudness = nrlz.*bandTonalLoudness;  % Equation 47 ECMA-418-2:2022

        % Section 6.2.8 Equation 48 ECMA-418-2:2022
        bandNoiseLoudness = acousticHMSNoiseRedLowPass_(bandLoudness, sampleRate1875) - bandTonalLoudness;  % specific loudness of non-tonal component in critical band
    
        % Store critical band results
        % ---------------------------
%         specificSNR(zBand, :, chan) = SNRlz1;  % specific time-dependent signal-noise-ratio in each critical band
%         specificLoudness(zBand, :, chan) = bandLoudness;  % specific time-dependent loudness of signal in each critical band
        specificTonalLoudness(zBand, :, chan) = bandTonalLoudness;  % specific time-dependent loudness of tonal component in each critical band
        specificNoiseLoudness(zBand, :, chan) = bandNoiseLoudness;   % specific time-dependent loudness of non-tonal component in each critical band
        specificTonalityFreqs(zBand, :, chan) = bandTonalFreqs;  % time-dependent frequency of tonal component in each critical band
        
    end

    % Calculation of specific tonality
    % --------------------------------
    % Section 6.2.8 Equation 49 ECMA-418-2:2022
    overallSNR = max(specificTonalLoudness, [], 1)./(1e-12 + sum(specificNoiseLoudness, 1));  % loudness signal-noise-ratio
    
    % Section 6.2.8 Equation 50 ECMA-418-2:2022
    crit = exp(-A*(overallSNR - B));
    ql = 1 - crit;  % sigmoidal scaling factor
    ql(crit >= 1) = 0;
    
    % Section 6.2.8 Equation 51 ECMA-418-2:2022#
    specificTonality = cal_T*cal_Tx*ql.*specificTonalLoudness;  % time-dependent specific tonality
    
    % Calculation of time-averaged specific tonality Section 6.2.9 ECMA-418-2:2022
    for zBand = 53:-1:1
        mask = specificTonality(zBand, :, chan) > 0.02;  % criterion Section 6.2.9 point 2
        mask(1:(58 - l_start)) = 0;  % criterion Section 6.2.9 point 1
        if l_start ~= 1
            mask(end + 1 - l_start:end) = 0;  % additional masking if 'ecma' is false
        end
        % Section 6.2.9 Equation 53 ECMA-418-2:2022
        specificTonalityAvg(zBand, 1, chan)...
            = sum(specificTonality(zBand, mask, chan), 2)./(nnz(mask) + 1e-12);
        specificTonalityAvgFreqs(zBand, 1, chan)...
            = sum(specificTonalityFreqs(zBand, mask, chan), 2)./(nnz(mask) + 1e-12);
    end

    % Calculation of total (non-specific) tonality Section 6.2.10
    % -----------------------------------------------------------
    % Further update can add the user input frequency range to determine
    % total tonality - not yet incorporated

    % Section 6.2.8 Equation 52 ECMA-418-2:2022
    % time (s) corresponding with results output
    t = (0:(size(specificTonality, 2) - 1))/sampleRate1875;

    % Section 6.2.10 Equation 61 ECMA-418-2:2022
    % Time-dependent total tonality
    [tonalityTimeVar(:, chan), zmax] = max(specificTonality(:, :, chan),...
                                           [], 1);
    for ll = size(specificTonalityFreqs, 2):-1:1
        tonalityTimeVarFreqs(ll, chan) = specificTonalityFreqs(zmax(ll), ll);
    end
    
    % Calculation of representative values Section 6.2.11 ECMA-418-2:2022
    % Time-averaged total tonality
    mask = tonalityTimeVar(:, chan) > 0.02;  % criterion Section 6.2.9 point 2
    mask(1:(58 - l_start)) = 0;    % criterion Section 6.2.9 point 1
    if l_start ~= 1
        mask(end + 1 - l_start:end) = 0;  % additional masking if 'ecma' is false
    end

    % Section 6.2.11 Equation 63 ECMA-418-2:2022
    % Time-averaged total tonality (note: epsilon is not applied here, according to the standard)
    tonalityAvg(chan) = sum(tonalityTimeVar(mask, chan))/nnz(mask);
    
    close(w)  % close waitbar

%% Output plotting

    % Plot figures
    % ------------
    if outplot
        cmap_plasma = load('cmap_plasma.txt');
        % Plot results
        chan_lab = chans(chan);
        fig = figure;
        tiledlayout(fig, 2, 1);
        movegui(fig, 'center');
        ax1 = nexttile(1);
        surf(ax1, t, bandCentreFreqs, specificTonality(:, :, chan),...
             'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        ax1.XLim = [t(1), t(end) + (t(2) - t(1))];
        ax1.YLim = [bandCentreFreqs(1), bandCentreFreqs(end)];
        ax1.CLim = [0, ceil(max(tonalityTimeVar(:, chan))*10)/10];
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
        set(get(h,'label'),'string', {'Specific Tonality,'; 'tu_{HMS}/Bark_{HMS}'});

        % Create A-weighting filter
        weightFilt = weightingFilter('A-weighting', sampleRate48k);
        % Filter signal to determine A-weighted time-averaged level
        pA = weightFilt(p_re(:, chan));
        LA = 20*log10(rms(pA)/2e-5);
        title(strcat(chan_lab,...
                     " signal sound pressure level =", {' '},...
                     num2str(round(LA,1)), "dB {\itL}_{Aeq}"),...
                     'FontWeight', 'normal', 'FontName', 'Arial');
        
        ax2 = nexttile(2);
        plot(ax2, t, tonalityTimeVar(:, chan), 'm', 'LineWidth', 1);
        ax2.XLim = [t(1), t(end) + (t(2) - t(1))];
        ax2.YLim = [0, ceil(max(tonalityTimeVar(:, chan))*10)/10];
        ax2.XLabel.String = "Time, s";
        ax2.YLabel.String = "Tonality, tu_{HMS}";
        ax2.XGrid = 'on';
        ax2.YGrid = 'on';
        ax2.FontName = 'Arial';
        ax2.FontSize = 12;
    end

end

% end of function