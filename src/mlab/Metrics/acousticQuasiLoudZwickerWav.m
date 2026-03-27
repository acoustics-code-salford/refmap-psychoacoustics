function loudness = acousticQuasiLoudZwickerWav(p, sampleRateIn, timeStep, axisN, soundField, adjustLoud, isoFilter, outPlot)
% loudness = acousticQuasiLoudZwickerWav(p, sampleRateIn, timeStep, axisN,
%                                        soundField, adjustLoud, isoFilter,
%                                        outPlot)
%
% Returns quasi-loudness using spectral elements of ISO 532-1 Zwicker
% loudness model for input pressure time-series. The temporal processing of
% the Zwicker model is disregarded. The input signal must be calibrated
% acoustic pressure (Pa).
%
% Optional modifications available comprise adjustments to the spectral
% levels to mimic the application of the ISO 532-3:2023 or
% ECMA-418-2:2025 outer-middle ear filter responses (in 1/3-octaves)
% instead of the ISO 532-1:2017 outer ear transmission, and
% approximated non-linear transformations more closely following
% ISO 532-3 or ECMA-418-2 instead of the ISO 532-1 original.
%
% Inputs
% ------
% p : vector or 2D matrix
%   Input signal as single mono or stereo audio (sound
%   pressure) signals.
%
% sampleRateIn : integer
%   Sample rate (frequency) of the input signal(s).
%
% timeStep : number
%   Time step value (seconds) used to calculate the time-dependent Leq
%   and loudness.
%
% axisN : integer (1 or 2, default: 1)
%   Time axis along which to calculate the loudness.
%
% soundField : keyword string (default: 'freeFrontal')
%   Determines whether the 'freeFrontal' or 'diffuse' field stages
%   are applied in the outer-middle ear filter, or 'noOuter'
%   omits this filtering stage.
%
% adjustLoud : keyword string (default: 'iso5321')
%   Indicates whether to apply adjustments for the outer-middle ear filter
%   response and loudness transformations from ISO 532-3:2023 ('iso5323')
%   or ECMA-418-2:2024 ('ecma4182').
%   The default option ('iso5321') applies no adjustment, so follows ISO
%   532-1 more closely (for closer agreement with Zwicker's model).
% 
% isoFilter : Boolean true/false (default: true)
%   Flag indicating whether to apply the ISO 532-1 1/3 octave band filters.
%   The ISO 532-1 filters are optimised for frequency response, but the
%   implementation is relatively slow. If false, a 6th order Butterworth 
%   1/3-octave filterbank is used instead, which is much faster (but less 
%   accurate, especially for low frequency bands).
%
% outPlot : Boolean true/false (default: false)
%   Flag indicating whether to generate a figure from the output.
%
% Returns
% -------
% loudness : structure
%   Contains the loudness output.
%
% loudness contains the following outputs:
%
% loudTDep :  vector or matrix
%   Time-dependent loudness arranged as [time(, channels)].
% 
% loudPowAvg : number or vector
%   Time-power-averaged loudness arranged as [loudness(, channels)].
%
% loud5pcEx : number or vector
%   95th percentile (5% exceeded) loudness
%   arranged as [loudness(, channels)].
%
% loudLevel :  vector or matrix
%   Time-dependent loudness level arranged as [time(, channels)].
%
% loudLvlPowAvg : number or vector
%   Time-power-averaged loudness level arranged as [loudness(, channels)].
%
% specLoudAvg : vector or matrix
%   Time-power-averaged specific loudness arranged as [bands(, channels)].
%
% specLoudAvg : matrix
%   Time-dependent specific loudness arranged as [time, bands(, channels)].
%
% barkAxis : vector
%   Critical band rates for specific loudness (0.1 dz intervals).
%
% freqInMid : vector
%   The 1/3-octave band exact mid-frequencies for the input spectra.
%
% freqInNom : vector
%   The 1/3-octave band nominal mid-frequencies for the input spectra.
%
% soundField : string
%   Identifies the soundfield type applied (the input argument soundField)
%
% timeOut : vector
%   Time (seconds) corresponding with time-dependent outputs.
%
% adjustLoud : string
%   Indicates which loudness model standard was adjusted for (=adjustLoud
%   input).
%
% Assumptions
% -----------
% The input is a 1/3-octave band unweighted sound level spectrum or
% series of sound level spectra in dB re 2e-5 Pa.
%
% References
% ----------
%
% Zwicker loudness is defined in ISO 532-1:2017. Modifications are based on
% ISO 532-3:2023 and ECMA-418-2:2025.
%
% Requirements
% ------------
% Signal Processing Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 23/04/2025
% Date last modified: 27/03/2026
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
        p (:, :) double {mustBeReal}
        sampleRateIn (1, 1) double {mustBePositive, mustBeInteger}
        timeStep (1, 1) double {mustBePositive}
        axisN (1, 1) {mustBeInteger, mustBeInRange(axisN, 1, 2)} = 1
        soundField (1, :) string {mustBeMember(soundField,...
                                                       {'freeFrontal',...
                                                        'diffuse',...
                                                        'noOuter'})} = 'freeFrontal'
        adjustLoud (1, :) string {mustBeMember(adjustLoud,...
                                               {'iso5321',...
                                                'iso5323',...
                                                'ecma4182'})} = 'iso5321'
        isoFilter (1, 1) {mustBeNumericOrLogical} = true
        outPlot (1, 1) {mustBeNumericOrLogical} = false
    end

%% Input checks

% Orientate signal
if axisN == 2
    p = p.';
end

% Ensure sampling rate is full range
resampledRate = 48e3;
if sampleRateIn ~= resampledRate  % Resample signal
    up = resampledRate/gcd(resampledRate, sampleRateIn);  % upsampling factor
    down = sampleRateIn/gcd(resampledRate, sampleRateIn);  % downsampling factor
    p_re = resample(p, up, down);  % apply resampling
else  % don't resample
    p_re = p;
end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Signal processing

% Get 1/3-octave filtered signals
if ~isoFilter
    % 6th order Butterworth filter bank
    octFiltBank = octaveFilterBank('Bandwidth', '1/3 octave',...
                                   'SampleRate', resampledRate,...
                                   'FrequencyRange', [25, 12600],...
                                   'FilterOrder', 6, ...
                                   'OctaveRatioBase', 10);

    signalFiltBank = octFiltBank(p_re);

else
    % ISO 532-1 compliant 6th order Chebyshev filter bank
    [signalFiltBank, ~] = iso5321_third_oct_filterbank(p_re, resampledRate, 1, [25, 12600]);
end

% Calculate Leq
[Leq, ~] = timeAveragedLevel(signalFiltBank, resampledRate, timeStep, 1);

loudness = acousticQuasiLoudZwicker(Leq, [25, 12600], timeStep, 2,...
                                    soundField, adjustLoud, outPlot);

% Recalculate time-aggregated outputs accounting for filter onset transient
if isoFilter
    % skip time
    timeSkip = 0.3;
    timeSkipIdx = ceil(timeSkip/timeStep) + 1;

    % power-averaged overall loudness
    loudness.loudPowAvg = power(mean(loudness.loudTDep(timeSkipIdx:end, :).^(1/log10(2)), 1), log10(2));

    % 95th percentile overall loudness
    loudness.loud5pcEx = prctile(loudness.loudTDep(timeSkipIdx:end, :), 95, 1);

    switch adjustLoud
        case 'iso5321'
            model = 'zwicker';
        case 'iso5323'
            model = 'mgs';
        case 'ecma4182'
            model = 'sottek';
    end
    loudness.loudLvlPowAvg = soneToPhon(loudness.loudPowAvg, model);
    loudness.loudLvl5pcEx = soneToPhon(loudness.loud5pcEx, model);

    % time-averaged specific loudness as a function of Bark number
    loudness.specLoudPowAvg = squeeze(power(mean(loudness.specLoud(timeSkipIdx:end, :, :).^(1/log10(2)), 1), log10(2)));

    % used for development purposes
    loudness.loudCorePowAvg = squeeze(power(mean(loudness.loudCore(timeSkipIdx:end, :, :).^(1/log10(2)), 1), log10(2)));
end

loudness.timeOut = linspace(0, size(loudness.loudTDep, 1)*timeStep...
                            - timeStep, size(loudness.loudTDep, 1)).';

end  % end of acousticQuasiLoudZwickerWav function