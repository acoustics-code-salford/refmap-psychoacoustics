function detectDiscount = acousticDetectDiscntFromLeqf(LeqTarget, LeqMasker, axisTarget, axisMasker, timeSkip, timeStep, freqMids, outPlot)
% detectDiscount = acousticDetectDiscntFromLeqf(LeqTarget, LeqMasker,
%                                               axisTarget, axisMasker,
%                                               timeSkip, timeSkip,
%                                               freqRange, outPlot)
%
% Returns detectability and discounted sound levels from input target
% source and masker signals based on the detectability model originally
% developed by Bolt, Beranek and Newman consulting engineers, and developed
% further by NASA, with discounting of target source sound levels using the
% technique developed by NASA.
%
% Fidell, S et al, 1974. Prediction of aural detectability of noise
% signals, Human Factors 16(4), 373-383.
% https://doi.org/10.1177/001872087401600405
%
% Sneddon, M et al, 2003. Laboratory study of the noticeability and
% annoyance of low signal-to-noise ratio sounds, Noise Control Engineering
% Journal, 51(5), 300-305.
% https://doi.org/10.3397/1.2839726
%
% Christian, A, 2021. A construction for the prediction of noise-induced
% annoyance in the presence of auditory masking, 181st Meeting of the
% Acoustical Society of America.
% https://ntrs.nasa.gov/citations/20210024824
%
% Rizzi, SA et al, 2024. Annoyance model assessments of urban air mobility
% vehicle operations. 30th AIAA/CEAS Aeroacoustics Conference, Rome, Italy,
% June 4-7, 2024. https://doi.org/10.2514/6.2024-3014
%
%
% Inputs
% ------
% signalTarget : vector or 2D matrix
%                the input target signal as single mono or stereo audio
%                (sound pressure) signals.
%
% signalMasker : vector or 2D matrix
%                the input masker signal(s) as single mono or stereo audio
%                (sound pressure) signals.
%
% axisTarget : integer (1 or 2, default: 1)
%              the time axis for the target signal(s) along which to determine
%              detection.
%
% axisMasker : integer (1 or 2, default: 1)
%              the time axis for the masker signal(s) along which to determine
%              detection.
%
% timeSkip : vector (default: [0, 0])
%            time (seconds) to skip from input signals for calculating
%            time-aggregated outputs. [startSkip, endSkip] ignores
%            starkSkip seconds of the start, and endSkip seconds of the
%            end.
%
% timeStep : number (default: 0.5)
%            the time window (seconds) to use for calculating target
%            detectability.
%
% freqMids : vector (default: [20, 20000])
%            the 1/3-octave band centre-frequencies for signal and
%            masker Leq inputs.
%
% outPlot : Boolean (default: false)
%           determines whether to plot outputs from the calculations
%
% Returns
% -------
% detectDiscount : structure
%                  contains the output
%
% detectDiscount contains the following outputs:
%
% dBSpecTarget : matrix
%                sound pressure level spectrogram for the input target
%                signal, with dimensions [timeOut, freqBands, targetChans]
%
% dBSpecMasker : matrix
%                sound pressure level spectrogram for the input masker
%                signal, with dimensions [timeOut, freqBands, maskerChans] 
%
% dBSpecDiscTarget : matrix
%                    sound pressure level spectrogram for the input target
%                    signal detectability-discounted, with dimensions
%                    [timeOut, freqBands, targetChans]
%
% dBATDepTarget : matrix or vector
%                 time-dependent A-weighted sound pressure level for the
%                 input target signal, with dimensions [timeOut, targetChans]
%
% dBATDepTarget : matrix or vector
%                 time-dependent A-weighted sound pressure level for the
%                 input masker signal, with dimensions [timeOut, targetChans]
%
% dBATDepDiscTarget : matrix or vector
%                     time-dependent A-weighted sound pressure level for the
%                     input target signal detectability-discounted, with
%                     dimensions [timeOut, targetChans]
%
% dBATDepDiscount : vector
%                   time-dependent detectability discount values for the
%                   A-weighted sound pressure levels of target signal vs
%                   masker signal, with dimensions [timeOut, targetChans]
%
% LAETarget : vector
%             overall A-weighted sound exposure level for each input target
%             signal channel
%
% LAEMasker : vector
%             overall A-weighted sound exposure level for each input masker
%             signal channel
%
% LAEDiscTarget : vector
%                 overall detectability-discounted A-weighted sound
%                 exposure level for each input target signal channel
%
% LAeqTarget : vector
%              overall A-weighted time-averaged sound level for each input
%              target signal channel
%
% LAeqMasker : vector
%              overall A-weighted time-averaged sound level for each input
%              masker signal channel
%
% LAeqDiscTarget : vector
%                  overall detectability-discounted A-weighted
%                  time-averaged sound level for each input target signal
%                  channel
%
% dBADiscount : vector
%               overall detectability discount for A-weighted levels for
%               each input target signal channel
%
% detectabilitydB : matrix
%                   detectability spectrogram for the input target signal
%                   vs masker signal (i.e., 10log_10[d']), with dimensions
%                   [timeOut, freqBands, targetChans]
%
% detectTDepMaxdB : matrix or vector
%                   band-maximum time-dependent detectability dB, with
%                   dimensions [timeOut, targetChans]
%
% detectTDepIntdB : matrix or vector
%                   band-integrated time-dependent detectability dB, with
%                   dimensions [timeOut, targetChans]
%
% detectMaxdB : vector
%               overall maximum detectability, dB, for each input target
%               signal channel
%
% detectIntdB : vector
%               overall intergated detectability, dB, for each input target
%               signal channel
%
% detectMaxPcdB : structure
%                 contains percent-based aggregated metrics for maximum
%                 detectability, comprising:
% Ex50 : vector
%        50% exceeded (median) dB, for each input target signal channel
%
% detectIntPcdB : structure
%                 contains percent-based aggregated metrics for maximum
%                 detectability, comprising:
% Ex50 : vector
%        50% exceeded (median) dB, for each input target signal channel
%
% freqBands : vector
%             1/3-octave band centre-frequencies for input freqBandRange
%
% timeOut : vector
%           the window centre times for the spectrograms
%
% Assumptions
% -----------
% The input signals are calibrated to units of acoustic pressure in Pascals
% (Pa).
%
% Requirements
% ------------
% Signal Processing Toolbox
% Audio Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 30/04/2025
% Date last modified: 30/04/2025
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
% Checked by:
% Date last checked:
%
%% Arguments validation;
    arguments (Input)
        LeqTarget (:, :) double {mustBeReal}
        LeqMasker (:, :) double {mustBeReal}
        axisTarget (1, 1) {mustBeInteger, mustBeInRange(axisTarget, 1, 2)} = 1
        axisMasker (1, 1) {mustBeInteger, mustBeInRange(axisMasker, 1, 2)} = 1
        timeSkip (1, 2) double {mustBeReal} = [0, 0]
        timeStep (1, 1) double {mustBePositive} = 0.5
        freqMids (1, 2) double {mustBeInRange(freqMids, 19, 20000)} = [20, 20000]
        outPlot {mustBeNumericOrLogical} = false
    end

%% Load path
addpath(genpath(fullfile("refmap-psychoacoustics", "src", "mlab")))

%% Input checks and resampling

% orient axes
if axisTarget == 2
    signalTarget = signalTarget.';
end

if axisMasker == 2
    signalMasker = signalMasker.';
end

% check input sizes and if necessary repeat masker
repMasker = 1;
targetChans = size(signalTarget, 2);
maskerChans = size(signalMasker, 2);
if targetChans ~= maskerChans
    % if masker is single channel, each target channel is assumed to be
    % matched with the masker (automatic broadcasting), otherwise, check if
    % target channels are a multiple of masker channels
    if maskerChans ~= 1
        if maskerChans > targetChans ||...
            mod(targetChans, maskerChans) ~= 0 
            error("Target and masker channel counts cannot be interpreted. Please ensure masker channel count is less than target and is divisible by target.")
        else
            % number of repeats for masker
            repMasker = targetChans/maskerChans;
        end
    end
end

% check input signal sampling frequencies are equal, otherwise resample to
% match higher rate
if sampleRateMasker > sampleRateTarget
    sampleRate = sampleRateMasker;
    up = sampleRate/gcd(sampleRate, sampleRateTarget);
    down = sampleRateTarget/gcd(sampleRate, sampleRateTarget);
    signalTarget = resample(signalTarget, up, down);
elseif sampleRateTarget > sampleRateMasker
    sampleRate = sampleRateTarget;
    up = sampleRate/gcd(sampleRate, sampleRateMasker);
    down = sampleRateMasker/gcd(sampleRate, sampleRateMasker);
    signalMasker = resample(signalMasker, up, down);
else
    sampleRate = sampleRateTarget;
end

% Check time skip
% Overall time (s)
T = size(signalTarget, 1)/sampleRate;

if sum(timeSkip) > T
    error("Error: total 'timeSkip' argument is larger than signal duration. Check inputs.")
end

% check frequency range
fl = min(freqRange);
fh = max(freqRange);

% ensure frequency range is suitable for signal sampling frequency
if fh > sampleRate/2.4
    fh = min(sampleRate/2.4, 20e3);
end

%% Define constants
b = 3;  % denominator for 1/b-octave definition

timeSteps = timeStep*sampleRate;  % timeStep in signal samples

efficiencyFactor= 0.3;  % \eta
targetDetect = 2;  % target d'
discountHalfPower = 3;  % \alpha, dB
detectKnee = 14;  % \delta, d' value at which 3 dB knee occurs
discountRate = 1;  % \rho, rate at which discount function dimishes with reducing detectability

% free-field frontal tone hearing thresholds for 18-25 year-olds with
% normal hearing from ISO 226:2023 (20 Hz - 12.5 kHz)
% hearThresholds226 = [78.1; 68.7; 59.5; 51.1; 44.0; 37.5; 31.5; 26.5; 22.1;...
%                      17.9; 14.4; 11.4; 8.6; 6.2; 4.4; 3.0; 2.2; 2.4; 3.5;...
%                      1.7; -1.3; -4.2; -6.0; -5.4; -1.5; 6.0; 12.6; 13.9; 12.3];

% free-field frontal tone hearing thresholds for 18-25 year-olds with
% normal hearing from ISO 389-7:2019 (16 - 18 kHz)
% hearThresholds3897 = [40.2; 70.4];

% free-field frontal tone hearing thresholds for 18-25 year-olds with
% normal hearing from ISO 226:2023 | ISO 389-7:2019 (20 - 18 kHz), with 20
% kHz band estimated from figure 1 in ISO 389-7:2019
% hearThresholds = [hearThresholds226; hearThresholds3897];

% diffuse field narrowband noise hearing thresholds for 18-25 year-olds with
% normal hearing from ISO 389-7:2019 (20 Hz - 16 kHz)
% NOTE: watch out for non-standard frequencies in the standard Table 1!
hearThresholdsDF3897 = [78.1; 68.7; 59.5; 51.1; 44.0; 37.5; 31.5; 26.5; 22.1;...
                        17.9; 14.4; 11.4; 8.4; 5.8; 3.8; 2.1; 1.0; 0.8; 1.9;...
                        0.5; -1.5; -3.1; -4.0; -3.8; -1.8; 2.5; 6.8; 8.4; 14.4;...
                        43.7];

% estimated full range diffuse field narrowband noise hearing thresholds
% for 18-25 year-olds with normal hearing 
hearThresholdsDF = [hearThresholdsDF3897; 70.4];

% 1/3-octave A-weighting dB values (20 Hz - 20 kHz)
Aweight = [-50.5; -44.7; -39.4; -34.6; -30.2; -26.2; -22.5; -19.1; -16.1;...
           -13.4; -10.9; -8.6; -6.6; -4.8; -3.2; -1.9; -0.8; 0.0; 0.6; 1.0;...
            1.2; 1.3; 1.2; 1.0; 0.5; -0.1; -1.1; -2.5; -4.3; -6.6; -9.3];


%% Signal processing

% convert timeSkip to timeStep indices
itimeSkip = timeSkip/timeStep;

% Calculate 1/3-octave power spectrograms
[spectroTarget, ~, ~] = poctave(signalTarget, sampleRate, 'spectrogram',...
                                'BandsPerOctave', b, 'WindowLength', timeSteps,...
                                'FrequencyLimits', [fl, fh]);

[spectroMasker, f, t] = poctave(signalMasker, sampleRate, 'spectrogram',...
                                'BandsPerOctave', b, 'WindowLength', timeSteps,...
                                'FrequencyLimits', [fl, fh]);

% Equivalent rectangular bandwidth for auditory filter (Glasberg & Moore,
% 1990, equation 3)
ERBandwidth = 24.7*(1 + 4.37/1000*f);

% Calculate 1/3 octave bandwidths according to IEC 61640-1:2014 / ANSI S1.11-2004
G10 = 10^(3/10);  % octave ratio coefficient (base-ten)
OctRatio = G10^(0.5/b);  % octave ratio

f1 = f/OctRatio;  % output range of exact lower band-edge frequencies
f2 = f*OctRatio;  % output range of exact upper band-edge frequencies
fBandWidth = f2 - f1;

% Calculate high frequency weighting "kick" factor
wKick = zeros(size(f));
wKick(f >= 2500) = -6*(log10(f( f>= 2500)/2500)).^2;

% Calculate detection efficiency
detectEfficiency = efficiencyFactor*sqrt(timeStep)*sqrt(fBandWidth).*sqrt(fBandWidth./ERBandwidth).*10.^(wKick/10);

% Calculate equivalent auditory system noise
ind = (-17:1:13).';  % range of frequency indices for hearing threshold bands
ind_vis = [ind; ind(end) + 1];  % used for visualisations
if mod(b, 1) == 0
    fm = G10.^(ind/b)*1000;
    fm_vis = G10.^(ind_vis/b)*1000;
else
    fm = G10.^((2*ind + 1)/(2*b))*1000;
    fm_vis = G10.^((2*ind_vis + 1)/(2*b))*1000;
end
fBandWidth_vis = fm_vis*OctRatio - fm_vis/OctRatio;

[~, il] = min(abs(fm - f(1)));  % find nearest lower exact frequency
[~, ih] = min(abs(fm - f(end)));  % find nearest higher exact frequency

eqAuditoryNoise = detectEfficiency.*(2e-5.*10.^(hearThresholdsDF(il:ih)/20)).^2/targetDetect;

% Calculate detectability and discount
detectability = detectEfficiency.*spectroTarget./(repmat(spectroMasker, 1, 1, repMasker) + eqAuditoryNoise);
detectabilitydB = 10*log10(detectability);
discountdB = discountHalfPower./(detectability./detectKnee).^discountRate;

% Calculate aggregated detectability
detectTDepMax = squeeze(max(detectability, [], 1));
detectTDepInt = squeeze(sqrt(sum(detectability.^2, 1)));
detectTDepMaxdB = 10*log10(detectTDepMax);
detectTDepIntdB = 10*log10(detectTDepInt);
detectMaxdB = 10*log10(max(detectTDepMax(1 + itimeSkip(1):end - itimeSkip(2), :), [], 1));
detectIntdB = 10*log10(timeStep*sum(detectTDepInt(1 + itimeSkip(1):end - itimeSkip(2), :), 1));
detectMaxPcdB.Ex50 = quantile(10*log10(detectTDepMax(1 + itimeSkip(1):end - itimeSkip(2), :)), 0.50, 1);
detectIntPcdB.Ex50 = quantile(10*log10(detectTDepInt(1 + itimeSkip(1):end - itimeSkip(2), :)), 0.50, 1);

% Calculate time-dependent spectra
dBSpecTarget = 20*log10(sqrt(spectroTarget)/2e-5);
dBSpecMasker = 20*log10(sqrt(spectroMasker)/2e-5);
dBSpecDiscTarget = dBSpecTarget - discountdB;

% A-weight time-dependent spectra
dBASpecTarget = dBSpecTarget + Aweight(il:ih);
dBASpecMasker = dBSpecMasker + Aweight(il:ih);
dBASpecDiscTarget = dBSpecDiscTarget + Aweight(il:ih);

% A-weighted time-dependent spectral power
powATarget = (2e-5*10.^(dBASpecTarget/20)).^2;
powAMasker = (2e-5*10.^(dBASpecMasker/20)).^2;
powADiscTarget = (2e-5*10.^(dBASpecDiscTarget/20)).^2;

% Aggregate spectral levels
dBATDepTarget = squeeze(20*log10(sqrt(sum(powATarget, 1))/2e-5));
dBATDepMasker = squeeze(20*log10(sqrt(sum(powAMasker, 1))/2e-5));
dBATDepDiscTarget = squeeze(20*log10(sqrt(sum(powADiscTarget, 1))/2e-5));

% Aggregate time-dependent levels
LAETarget = squeeze(20*log10(sqrt(timeStep*sum(sum(powATarget(1 + itimeSkip(1):end - itimeSkip(2), :, :), 1), 2))/2e-5));
LAEMasker = squeeze(20*log10(sqrt(timeStep*sum(sum(powAMasker(1 + itimeSkip(1):end - itimeSkip(2), :, :), 1), 2))/2e-5));
LAEDiscTarget = squeeze(20*log10(sqrt(timeStep*sum(sum(powADiscTarget(1 + itimeSkip(1):end - itimeSkip(2), :, :), 1), 2))/2e-5));
LAeqTarget = LAETarget - 10*log10(T);
LAeqMasker = LAEMasker - 10*log10(T);
LAeqDiscTarget = LAEDiscTarget - 10*log10(T);

% Time-dependent and overall dBA discounts
dBATDepDiscount = dBATDepTarget - dBATDepDiscTarget;
dBADiscount = LAETarget - LAEDiscTarget;

%% Assign outputs
detectDiscount.dBSpecTarget = dBSpecTarget;
detectDiscount.dBSpecMasker = dBSpecMasker;
detectDiscount.dBSpecDiscTarget = dBSpecDiscTarget;
detectDiscount.dBATDepTarget = dBATDepTarget;
detectDiscount.dBATDepMasker = dBATDepMasker;
detectDiscount.dBATDepDiscTarget = dBATDepDiscTarget;
detectDiscount.dBATDepDiscount = dBATDepDiscount;
detectDiscount.LAETarget = LAETarget;
detectDiscount.LAEMasker = LAEMasker;
detectDiscount.LAEDiscTarget = LAEDiscTarget;
detectDiscount.LAeqTarget = LAeqTarget;
detectDiscount.LAeqMasker = LAeqMasker;
detectDiscount.LAeqDiscTarget = LAeqDiscTarget;
detectDiscount.dBADiscount = dBADiscount;
detectDiscount.detectabilitydB = detectabilitydB;
detectDiscount.detectTDepMaxdB = detectTDepMaxdB;
detectDiscount.detectTDepIntdB = detectTDepIntdB;
detectDiscount.detectMaxdB = detectMaxdB;
detectDiscount.detectIntdB = detectIntdB;
detectDiscount.detectMaxPcdB.Ex50 = detectMaxPcdB.Ex50;
detectDiscount.detectIntPcdB.Ex50 = detectIntPcdB.Ex50;
detectDiscount.freqBands = f;
detectDiscount.timeOut = t.';

%% Plotting

if outPlot
    % Plot figures
    % ------------
    cmap_viridis = load('cmap_viridis.txt');
    cmap_magma = load('cmap_magma.txt');
    
    % number of channels for plotting
    if repMasker == 1
        chanMap = [1:targetChans; ones(1, targetChans)];
    else
        chanMap = [1:targetChans; repmat(1:maskerChans, 1, repMasker)];
    end

    % loop through channels with tiled figure for each
    for iiPlot = 1:targetChans
        targChan = chanMap(1, iiPlot);
        maskChan = chanMap(2, iiPlot);

        fig = figure;
        set(fig, 'position', [300, 200, 1300, 600])
        tiledlayout(fig, 2, 3);
        movegui(fig, 'center');

        ax1 = nexttile(1);
        surf(ax1, [t, t(end) + timeStep] - timeStep/2, fm_vis - fBandWidth_vis/2,...
             arrangeSpectro(dBSpecTarget(:, :, targChan)),...
             EdgeColor='none', FaceColor='interp');
        view(2);
        set(ax1, 'XLim', [t(1) - timeStep/2, t(end) - timeStep/2 + timeStep],...
            'YScale', 'log', 'YLim', [f1(1), f2(end)],...
            'YTick', [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],...
            'YTickLabel', ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]);
        ax1.YLabel.String = "Frequency, Hz";
        ax1.XLabel.String = "Time, s";
        colormap(ax1, cmap_viridis);
        cbar = colorbar; clim([0, max(dBSpecTarget(:, :, targChan), [], 'all')]);
        cbar.Label.String = "dB re 2e-5 Pa";
        ax1.Title.String = "Target spectrogram";
        ax1.TitleFontWeight = "normal";

        ax2 = nexttile(2);
        surf(ax2, [t, t(end) + timeStep] - timeStep/2, fm_vis - fBandWidth_vis/2,...
             arrangeSpectro(detectabilitydB(:, :, targChan)),...
             EdgeColor='none', FaceColor='interp');
        view(2);
        set(ax2, 'XLim', [t(1) - timeStep/2, t(end) - timeStep/2 + timeStep],...
            'YScale', 'log', 'YLim', [f1(1), f2(end)],...
            'YTick', [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],...
            'YTickLabel', ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]);
        ax2.YLabel.String = "Frequency, Hz";
        ax2.XLabel.String = "Time, s";
        colormap(ax2, cmap_magma);
        cbar = colorbar; clim([0, max(detectabilitydB(:, :, targChan), [], 'all')]);
        cbar.Label.String = "Detectability 10log_{10}{\it d'}, dB";
        ax2.Title.String = "Masked target detectability";
        ax2.TitleFontWeight = "normal";

        ax3 = nexttile(3);
        plot(ax3, t, dBATDepTarget(:, targChan), 'color', cmap_magma(34, :), 'DisplayName', "Target")
        hold on
        plot(ax3, t, dBATDepMasker(:, maskChan), 'color', cmap_magma(166, :), 'DisplayName', "Masker")
        plot(ax3, t, dBATDepDiscTarget(:, targChan), 'color', cmap_magma(100, :),...
            'LineStyle', ':', 'LineWidth', 2, 'DisplayName', "Target discounted")
        hold off
        set(ax3, 'XLim', [t(1) - timeStep/2, t(end) - timeStep/2 + timeStep],...
            'YLim', [min(dBATDepMasker(:, maskChan), [], 'all') - 10,...
                     max(max(dBATDepTarget(:, targChan), [], 'all'),...
                     max(dBATDepMasker(:, maskChan), [], 'all'))*1.05],...
            'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.075,...
            'GridLineStyle', '--', 'GridLineWidth', 0.25);
        ax3.YLabel.String = "dB(A) re 2e-5 Pa";
        ax3.XLabel.String = "Time, s";
        ax3.Title.String = "Time-dependent levels";
        ax3.TitleFontWeight = "normal";
        legend(ax3, 'Location', 'eastoutside');
        
        ax4 = nexttile(4);
        surf(ax4, [t, t(end) + timeStep] - timeStep/2, fm_vis - fBandWidth_vis/2,...
             arrangeSpectro(dBSpecMasker(:, :, targChan)),...
             EdgeColor='none', FaceColor='interp');
        view(2);
        set(ax4, 'XLim', [t(1) - timeStep/2, t(end) - timeStep/2 + timeStep],...
            'YScale', 'log', 'YLim', [f1(1), f2(end)],...
            'YTick', [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],...
            'YTickLabel', ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]);
        ax4.YLabel.String = "Frequency, Hz";
        ax4.XLabel.String = "Time, s";
        colormap(ax4, cmap_viridis);
        cbar = colorbar; clim([0, max(dBSpecMasker(:, :, targChan), [], 'all')]);
        cbar.Label.String = "dB re 2e-5 Pa";
        ax4.Title.String = "Masker spectrogram";
        ax4.TitleFontWeight = "normal";

        ax5 = nexttile(5);
        surf(ax5, [t, t(end) + timeStep] - timeStep/2, fm_vis - fBandWidth_vis/2,...
             arrangeSpectro(dBSpecDiscTarget(:, :, targChan)),...
             EdgeColor='none', FaceColor='interp');
        view(2);
        set(gca, 'XLim', [t(1) - timeStep/2, t(end) - timeStep/2 + timeStep],...
            'YScale', 'log', 'YLim', [f1(1), f2(end)],...
            'YTick', [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],...
            'YTickLabel', ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]);
        ax5.YLabel.String = "Frequency, Hz";
        ax5.XLabel.String = "Time, s";
        colormap(ax5, cmap_magma);
        cbar = colorbar; clim([0, max(dBSpecDiscTarget(:, :, targChan), [], 'all')]);
        cbar.Label.String = "Detectability-discounted \newline     level, dB re 2e-5 Pa";
        ax5.Title.String = "Target detectability-discounted spectrogram";
        ax5.TitleFontWeight = "normal";

        ax6 = nexttile(6);
        levelVals = [LAETarget(targChan), LAEMasker(maskChan), LAEDiscTarget(targChan)];
        labelVals = num2cell(round(levelVals, 1));
        labelCats = ["Target", "Masker", "Discount. targ."];
        % a trick using stacked bar to get the legend mapped
        b = bar(ax6,...
                diag(levelVals, 0),...
                'stacked', 'FaceColor', 'flat');
        set(b, {'FaceColor'}, {cmap_magma(34, :); cmap_magma(166, :); cmap_magma(100, :)})
        set(ax6, 'YLim', [min(levelVals, [], 'all')/1.25, max(levelVals)*1.1])
        ax6.YLabel.String = "{\it L}_{AE}, dB re 2e-5 Pa";
        ax6.XTickLabel = [];
        lg6 = legend(ax6, labelCats, 'Location','eastoutside');
        lg6.Direction = 'normal';
        % data labels
        y = sum(reshape(cell2mat(get(b', 'YData')),size(b, 2), []), 1); 
        x = unique(cell2mat(get(b', 'XData')),'stable');
        offset = range(ylim)*.1; 
        text(x, y - offset, labelVals, 'HorizontalAlignment', 'Center',...
             'VerticalAlignment', 'bottom', 'Color', 'w');
        ax6.Title.String = "Overall levels";
        ax6.TitleFontWeight = "normal";
    end
end % end of if branch for plotting outputs

function returnXft = arrangeSpectro(targetXft)
% returnXft = arrangeSpectro(targetXft)
%
% Return 2D matrix arranged for plotting spectrogram-like plots with fully
% rendered cells at boundaries. Input and output matrices have shape
% [freqs, time]

returnXft = [targetXft(:, :), targetXft(:, end);
             targetXft(end, :), targetXft(end, end)];

% end of function