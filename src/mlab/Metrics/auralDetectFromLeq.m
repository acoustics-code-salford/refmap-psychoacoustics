function detectability = auralDetectFromLeq(leqtSpecTarget, leqtSpecMasker, timeStep, timeSkip, nOct, freqRange, outPlot)
% detectability = auralDetectFromLeq(leqtSpecTarget, leqtSpecMasker, timeStep,
%                                     timeSkip, nOct, freqRange, outPlot)
%
% Returns aural detectability and discounted sound levels from input target
% source and masker spectral sound level time series [LZeq(t)] data, based
% on the detectability model originally developed by Bolt, Beranek and 
% Newman consulting engineers, and developed further by NASA, with 
% discounting of target source sound levels using the technique developed 
% by NASA (see References).
%
% Inputs
% ------
% leqtSpecTarget : 2D or 3D matrix
%   Target signal LZeq(t, f) spectrogram [time, freq_bands(, channels)].
% 
% leqtSpecMasker : 2D or 3D matrix
%   Masker signal LZeq(t, f) spectrogram [time, freq_bands(, channels)].
%
% timeStep : number (default: 0.5)
%   Time window (seconds) to use for calculating target detectability.
%
% timeSkip : vector (default: [0, 0])
%   Time (seconds) to skip from input signals for calculating
%   time-aggregated outputs. [startSkip, endSkip] ignores
%   starkSkip seconds of the start, and endSkip seconds of the end.
%
% nOct : integer (default: 3)
%   Number of fractional-octave bands to use (e.g., default 3 = 1/3-octave).
%   Note that fractional-octave bands other than 1/3-octave have not yet
%   been implemented.
%
% freqRange : vector (default: [19, 20000])
%   Frequency range over which to determine detection and discounted 
%   spectra (1/3-octave band centre-frequencies within this range will be 
%   included)
%
% outPlot : Boolean (default: false)
%   Determines whether to plot outputs from the calculations.
% 
% Returns
% -------
% detectability : structure
%   Contains the output.
%
% detectability contains the following outputs:
%
% leqtSpecTarget : 2D or 3D matrix
%   Target signal LZeq(t, f) spectrogram [time, freq_bands(, channels)].
% 
% leqtSpecMasker : 2D or 3D matrix
%   Masker signal LZeq(t, f) spectrogram [time, freq_bands(, channels)].
%
% leqtSpecDiscTarget : 2D or 3D matrix
%   Detectability-discounted LZeq(t, f) spectrogram for the input target
%   signal, with dimensions [timeOut, freqBands, targetChans].
%
% lAeqtTarget : matrix or vector
%   Time-dependent A-weighted LAeq(t) for the input target signal, 
%   with dimensions [timeOut, targetChans].
%
% lAeqtMasker : matrix or vector
%   Time-dependent A-weighted LAeq(t) for the input masker signal,
%   with dimensions [timeOut, targetChans].
%
% lAeqtDiscTarget : matrix or vector
%   Time-dependent detectability-discounted A-weighted sound pressure level
%   for the input target signal, with dimensions [timeOut, targetChans].
%
% lAeqTarget : vector
%   Overall A-weighted time-averaged sound level for each input
%   target signal channel.
%
% lAeqMasker : vector
%   Overall A-weighted time-averaged sound level for each input
%   masker signal channel.
%
% lAeqDiscTarget : vector
%   Overall detectability-discounted A-weighted time-averaged sound level 
%   for each input target signal channel.
%
% lAETarget : vector
%   Overall A-weighted sound exposure level for each input target signal 
%   channel.
%
% lAEMasker : vector
%   Overall A-weighted sound exposure level for each input masker signal 
%   channel.
%
% lAEDiscTarget : vector
%   Overall detectability-discounted A-weighted sound exposure level for 
%   each input target signal channel.
%
% dBADiscount : vector
%   Overall detectability discount for A-weighted levels for each input 
%   target signal channel.
%
% detectabilitydB : matrix
%   Detectability spectrogram for the input target signal
%   vs masker signal (i.e., 10log_10[d']), with dimensions
%   [timeOut, freqBands, targetChans].
%
% detectTMaxdB : matrix or vector
%   Band-maximum time-dependent detectability dB, with dimensions 
%   [timeOut, targetChans].
%
% detectTIntdB : matrix or vector
%   Band-integrated time-dependent detectability dB, with dimensions 
%   [timeOut, targetChans].
%
% detectMaxdB : vector
%   Overall maximum detectability, dB, for each input target signal
%   channel.
%
% detectIntdB : vector
%   Overall integrated detectability, dB, for each input target signal 
%   channel.
% 
% detectMaxIntdB : vector
%   Overall spectral-maximum, time-integrated detectability, dB, for each 
%   input target signal channel.
% 
% detectIntMaxdB : vector
%   Overall spectral-integrated, time-maximum detectability, dB, for each 
%   input target signal channel.
%
% detectMaxPcdB : structure
%   Contains percent-based time-aggregated metrics for
%   spectral-maximum detectability, comprising:
%
% Ex50 : vector
%   50% exceeded (median) dB, for each input target signal channel.
%
% detectIntPcdB : structure
%   Contains percent-based time-aggregated metrics for
%   spectral-integrated detectability, comprising:
%
% Ex50 : vector
%   50% exceeded (median) dB, for each input target signal channel.
%
% freqBands : vector
%   1/3-octave band centre-frequencies for input freqRange.
%
% timeOut : vector
%   Window start times for the time-dependent outputs.
%
% Assumptions
% -----------
% The input spectra are orientated with time along axis 1 and frequency 
% bands on axis 2.
%
% References
% ----------
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
% Date created: 05/11/2024
% Date last modified: 23/03/2026
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
        leqtSpecTarget (:, :, :) double {mustBeReal}
        leqtSpecMasker (:, :, :) double {mustBeReal}
        timeStep (1, 1) double {mustBePositive} = 0.5
        timeSkip (1, 2) double {mustBeReal, mustBeNonnegative} = [0, 0]
        nOct (1, 1) {mustBeInteger, mustBePositive} = 3
        freqRange (1, 2) double {mustBeInRange(freqRange, 19, 20000)} = [19, 20000]
        outPlot {mustBeNumericOrLogical} = false
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Input checks and resampling
% check input lengths match
if size(leqtSpecTarget, 1) ~= size(leqtSpecMasker, 1) || size(leqtSpecTarget, 2) ~= size(leqtSpecMasker, 2)
    error("Error: The lengths of the input spectra on the time or frequency band axes do not match.")
end

% check input channel sizes and if necessary repeat masker
repMasker = 1;
if ndims(leqtSpecTarget) > 2 %#ok<ISMAT>
    targetChans = size(leqtSpecTarget, 3);
else
    targetChans = size(leqtSpecTarget, 2);
end
if ndims(leqtSpecMasker) > 2 %#ok<ISMAT>
    maskerChans = size(leqtSpecMasker, 3);
else
    maskerChans = size(leqtSpecMasker, 2);
end

if targetChans ~= maskerChans
    % if masker is single channel, each target channel is assumed to be
    % matched with the masker (automatic broadcasting), otherwise, check if
    % target channels are a multiple of masker channels
    if maskerChans ~= 1
        if maskerChans > targetChans ||...
            mod(targetChans, maskerChans) ~= 0 
            error("Error: Target and masker channel counts cannot be interpreted. Please ensure masker channel count is less than target and is divisible by target.")
        else
            % number of repeats for masker
            repMasker = targetChans/maskerChans;
        end
    end
end

% Check time skip
% Overall time (s)
T = size(leqtSpecTarget, 1)*timeStep;
if sum(timeSkip) > T
    error("Error: total 'timeSkip' argument is larger than signal duration. Check inputs.")
end

%% Define constants
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
hearThresholdsDF3897 = [78.1, 68.7, 59.5, 51.1, 44.0, 37.5, 31.5, 26.5, 22.1,...
                        17.9, 14.4, 11.4, 8.4, 5.8, 3.8, 2.1, 1.0, 0.8, 1.9,...
                        0.5, -1.5, -3.1, -4.0, -3.8, -1.8, 2.5, 6.8, 8.4, 14.4,...
                        43.7];

% estimated full range diffuse field narrowband noise hearing thresholds
% for 18-25 year-olds with normal hearing 
hearThresholdsDF = [hearThresholdsDF3897, 70.4];

% 1/3-octave A-weighting dB values (20 Hz - 20 kHz)
Aweight = [-50.5, -44.7, -39.4, -34.6, -30.2, -26.2, -22.5, -19.1, -16.1,...
           -13.4, -10.9, -8.6, -6.6, -4.8, -3.2, -1.9, -0.8, 0.0, 0.6, 1.0,...
            1.2, 1.3, 1.2, 1.0, 0.5, -0.1, -1.1, -2.5, -4.3, -6.6, -9.3];

% Calculate 1/3 octave bandwidths according to IEC 61640-1:2014 / ANSI S1.11-2004
g10 = 10^(3/10);  % octave ratio coefficient (base-ten)
octRatio = g10^(0.5/nOct);  % octave ratio
fref = 1e3;  % reference frequency
indWide = -15*nOct:1:15*nOct;  % large range of indices to calculate frequencies

if mod(nOct, 1) == 0
    f = g10.^(indWide./nOct)*fref;
else
    f = g10.^((2*indWide + 1)./(2*nOct))*fref;
end

% find range of relevance
[~, ilow] = min(abs(f - min(freqRange)));  % index of highest f in frange
[~, ihigh] = min(abs(f - max(freqRange)));  % index of lowest f in frange

% band mid frequencies
fmid = f(ilow:ihigh);

if size(leqtSpecTarget, 2) ~= length(fmid) || size(leqtSpecMasker, 2) ~= length(fmid)
    error("Error: the input range of frequency bands (" + num2str(length(fmid)) +...
          ") does not match the length of the input signal frequency band axes (target: "...
          + num2str(size(leqtSpecTarget, 2)) + ", masker: " + num2str(size(leqtSpecMasker, 2)) + ").")
end

f1 = fmid/octRatio;  % output range of exact lower band-edge frequencies
f2 = fmid*octRatio;  % output range of exact upper band-edge frequencies
fBandWidth = f2 - f1;

% convert timeSkip to timeStep indices
itimeSkip = timeSkip/timeStep;

% Equivalent rectangular bandwidth for auditory filter (Glasberg & Moore,
% 1990, equation 3)
ERBandwidth = 24.7*(1 + 4.37/1000*fmid);

%% Signal processing

% calculate signal power spectral time series from Leqt
powTarget = (2e-5*10.^(leqtSpecTarget/20)).^2;
powMasker = (2e-5*10.^(leqtSpecMasker/20)).^2;

% Calculate high frequency weighting "kick" factor
wKick = zeros(size(fmid));
wKick(fmid >= 2500) = -6*(log10(fmid(fmid >= 2500)/2500)).^2;

% Calculate detection efficiency
detectEfficiency = efficiencyFactor*sqrt(timeStep)*sqrt(fBandWidth).*sqrt(fBandWidth./ERBandwidth).*10.^(wKick/10);

% Calculate equivalent auditory system noise
ind = (-17:1:13).';  % range of frequency indices for hearing threshold bands
ind_vis = [ind; ind(end) + 1];  % used for visualisations
if mod(nOct, 1) == 0
    fm = g10.^(ind/nOct)*1000;
    fm_vis = g10.^(ind_vis/nOct)*1000;
else
    fm = g10.^((2*ind + 1)/(2*nOct))*1000;
    fm_vis = g10.^((2*ind_vis + 1)/(2*nOct))*1000;
end
fBandWidth_vis = fm_vis*octRatio - fm_vis/octRatio;

[~, il] = min(abs(fm - fmid(1)));  % find nearest lower exact frequency
[~, ih] = min(abs(fm - fmid(end)));  % find nearest higher exact frequency

eqAuditoryNoise = detectEfficiency.*(2e-5.*10.^(hearThresholdsDF(il:ih)/20)).^2/targetDetect;

% Calculate detectability and discount
detectSpecT = detectEfficiency.*powTarget./(repmat(powMasker, 1, 1, repMasker) + eqAuditoryNoise);
detectSpecTdB = 10*log10(detectSpecT);
discountdB = discountHalfPower./(detectSpecT./detectKnee).^discountRate;
leqtSpecDiscTarget = leqtSpecTarget - discountdB;

% A-weight time-dependent spectra
lAeqSpecTarget = leqtSpecTarget + Aweight(il:ih);
lAeqSpecMasker = leqtSpecMasker + Aweight(il:ih);
lAeqSpecDiscTarget = leqtSpecDiscTarget + Aweight(il:ih);

% A-weighted time-dependent spectral power
powATarget = (2e-5*10.^(lAeqSpecTarget/20)).^2;
powAMasker = (2e-5*10.^(lAeqSpecMasker/20)).^2;
powADiscTarget = (2e-5*10.^(lAeqSpecDiscTarget/20)).^2;

% Calculate spectral aggregation of time-dependent values
detectTMax = squeeze(max(detectSpecT, [], 1));
detectTInt = squeeze(sqrt(sum(detectSpecT.^2, 1)));
lAeqtTarget = squeeze(20*log10(sqrt(sum(powATarget, 2))/2e-5));
lAeqtMasker = squeeze(20*log10(sqrt(sum(powAMasker, 2))/2e-5));
lAeqtDiscTarget = squeeze(20*log10(sqrt(sum(powADiscTarget, 2))/2e-5));
% dBAtDiscount = lAeqtTarget - lAeqtDiscTarget;
detectTMaxdB = 10*log10(detectTMax);
detectTIntdB = 10*log10(detectTInt);

% Calculate time-aggregated values
detectMaxdB = 10*log10(max(detectTMax(1 + itimeSkip(1):end - itimeSkip(2), :), [], 1));
detectIntdB = 10*log10(timeStep*sum(detectTInt(1 + itimeSkip(1):end - itimeSkip(2), :), 1));
detectMaxIntdB = 10*log10(timeStep*sum(detectTMax(1 + itimeSkip(1):end - itimeSkip(2), :), 1));
detectIntMaxdB = 10*log10(max(detectTInt(1 + itimeSkip(1):end - itimeSkip(2), :), [], 1));
detectMaxPcdB.Ex50 = quantile(10*log10(detectTMax(1 + itimeSkip(1):end - itimeSkip(2), :)), 0.50, 1);
detectIntPcdB.Ex50 = quantile(10*log10(detectTInt(1 + itimeSkip(1):end - itimeSkip(2), :)), 0.50, 1);
lAETarget = squeeze(20*log10(sqrt(timeStep*sum(sum(powATarget(1 + itimeSkip(1):end - itimeSkip(2), :, :), 2), 1))/2e-5));
lAEMasker = squeeze(20*log10(sqrt(timeStep*sum(sum(powAMasker(1 + itimeSkip(1):end - itimeSkip(2), :, :), 2), 1))/2e-5));
lAEDiscTarget = squeeze(20*log10(sqrt(timeStep*sum(sum(powADiscTarget(1 + itimeSkip(1):end - itimeSkip(2), :, :), 2), 1))/2e-5));
lAeqTarget = lAETarget - 10*log10(T);
lAeqMasker = lAEMasker - 10*log10(T);
lAeqDiscTarget = lAEDiscTarget - 10*log10(T);
dBADiscount = lAETarget - lAEDiscTarget;

% Output time
timeOut = timeStep*(0:size(leqtSpecTarget, 1) - 1).';

%% Assign outputs
detectability.leqtSpecTarget = leqtSpecTarget;
detectability.leqtSpecMasker = leqtSpecMasker;
detectability.leqtSpecDiscTarget = leqtSpecDiscTarget;
detectability.lAeqtTarget = lAeqtTarget;
detectability.lAeqtMasker = lAeqtMasker;
detectability.lAeqtDiscTarget = lAeqtDiscTarget;
detectability.lAeqTarget = lAeqTarget;
detectability.lAeqMasker = lAeqMasker;
detectability.lAeqDiscTarget = lAeqDiscTarget;
detectability.lAETarget = lAETarget;
detectability.lAEMasker = lAEMasker;
detectability.lAEDiscTarget = lAEDiscTarget;
detectability.dBADiscount = dBADiscount;
detectability.detectabilitydB = detectSpecTdB;
detectability.detectTMaxdB = detectTMaxdB;
detectability.detectTIntdB = detectTIntdB;
detectability.detectMaxdB = detectMaxdB;
detectability.detectIntdB = detectIntdB;
detectability.detectMaxIntdB = detectMaxIntdB;
detectability.detectIntMaxdB = detectIntMaxdB;
detectability.detectMaxPcdB.Ex50 = detectMaxPcdB.Ex50;
detectability.detectIntPcdB.Ex50 = detectIntPcdB.Ex50;
detectability.freqBands = fmid;
detectability.timeOut = timeOut;

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
        surf(ax1, [timeOut; timeOut(end) + timeStep], fm_vis - fBandWidth_vis/2,...
             arrangeSpectro(leqtSpecTarget(:, :, targChan)).',...
             EdgeColor='none', FaceColor='interp');
        view(2);
        set(ax1, 'XLim', [timeOut(1); timeOut(end) + timeStep],...
            'YScale', 'log', 'YLim', [f1(1), f2(end)],...
            'YTick', [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],...
            'YTickLabel', ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]);
        ax1.YLabel.String = "Frequency, Hz";
        ax1.XLabel.String = "Time, s";
        colormap(ax1, cmap_viridis);
        cbar = colorbar; clim([0, max(leqtSpecTarget(:, :, targChan), [], 'all')]);
        cbar.Label.String = "dB re 2e-5 Pa";
        ax1.Title.String = "Target {\it L}_{Zeq} spectrogram";
        ax1.TitleFontWeight = "normal";

        ax2 = nexttile(2);
        surf(ax2, [timeOut; timeOut(end) + timeStep], fm_vis - fBandWidth_vis/2,...
             arrangeSpectro(detectSpecTdB(:, :, targChan)).',...
             EdgeColor='none', FaceColor='interp');
        view(2);
        set(ax2, 'XLim', [timeOut(1); timeOut(end) + timeStep],...
            'YScale', 'log', 'YLim', [f1(1), f2(end)],...
            'YTick', [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],...
            'YTickLabel', ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]);
        ax2.YLabel.String = "Frequency, Hz";
        ax2.XLabel.String = "Time, s";
        colormap(ax2, cmap_magma);
        cbar = colorbar; clim([0, max(detectSpecTdB(:, :, targChan), [], 'all')]);
        cbar.Label.String = "Detectability 10log_{10}{\it d'}, dB";
        ax2.Title.String = "Masked target detectability";
        ax2.TitleFontWeight = "normal";

        ax3 = nexttile(3);
        plot(ax3, timeOut, lAeqtTarget(:, targChan), 'color', cmap_magma(34, :), 'DisplayName', "Target")
        hold on
        plot(ax3, timeOut, lAeqtMasker(:, maskChan), 'color', cmap_magma(166, :), 'DisplayName', "Masker")
        plot(ax3, timeOut, lAeqtDiscTarget(:, targChan), 'color', cmap_magma(100, :),...
            'LineStyle', ':', 'LineWidth', 2, 'DisplayName', "Target discounted")
        hold off
        set(ax3, 'XLim', [timeOut(1); timeOut(end) + timeStep],...
            'YLim', [min(lAeqtMasker(:, maskChan), [], 'all') - 10,...
                     max(max(lAeqtTarget(:, targChan), [], 'all'),...
                     max(lAeqtMasker(:, maskChan), [], 'all'))*1.05],...
            'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.075,...
            'GridLineStyle', '--', 'GridLineWidth', 0.25);
        ax3.YLabel.String = "dB(A) re 2e-5 Pa";
        ax3.XLabel.String = "Time, s";
        ax3.Title.String = "Time-dependent levels";
        ax3.TitleFontWeight = "normal";
        legend(ax3, 'Location', 'eastoutside');
        
        ax4 = nexttile(4);
        surf(ax4, [timeOut; timeOut(end) + timeStep], fm_vis - fBandWidth_vis/2,...
             arrangeSpectro(leqtSpecMasker(:, :, targChan)).',...
             EdgeColor='none', FaceColor='interp');
        view(2);
        set(ax4, 'XLim', [timeOut(1); timeOut(end) + timeStep],...
            'YScale', 'log', 'YLim', [f1(1), f2(end)],...
            'YTick', [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],...
            'YTickLabel', ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]);
        ax4.YLabel.String = "Frequency, Hz";
        ax4.XLabel.String = "Time, s";
        colormap(ax4, cmap_viridis);
        cbar = colorbar; clim([0, max(leqtSpecMasker(:, :, targChan), [], 'all')]);
        cbar.Label.String = "dB re 2e-5 Pa";
        ax4.Title.String = "Masker {\it L}_{Zeq} spectrogram";
        ax4.TitleFontWeight = "normal";

        ax5 = nexttile(5);
        surf(ax5, [timeOut; timeOut(end) + timeStep], fm_vis - fBandWidth_vis/2,...
             arrangeSpectro(leqtSpecDiscTarget(:, :, targChan)).',...
             EdgeColor='none', FaceColor='interp');
        view(2);
        set(gca, 'XLim', [timeOut(1); timeOut(end) + timeStep],...
            'YScale', 'log', 'YLim', [f1(1), f2(end)],...
            'YTick', [31.5, 63, 125, 250, 500, 1e3, 2e3, 4e3, 8e3, 16e3],...
            'YTickLabel', ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]);
        ax5.YLabel.String = "Frequency, Hz";
        ax5.XLabel.String = "Time, s";
        colormap(ax5, cmap_magma);
        cbar = colorbar; clim([0, max(leqtSpecDiscTarget(:, :, targChan), [], 'all')]);
        cbar.Label.String = "Detectability-discounted \newline     level, dB re 2e-5 Pa";
        ax5.Title.String = "Target detectability-discounted {\it L}_{Zeq} spectrogram";
        ax5.TitleFontWeight = "normal";

        ax6 = nexttile(6);
        levelVals = [lAETarget(targChan), lAEMasker(maskChan), lAEDiscTarget(targChan)];
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
% rendered cells at boundaries.

returnXft = [targetXft(:, :), targetXft(:, end);
             targetXft(end, :), targetXft(end, end)];

% end of function