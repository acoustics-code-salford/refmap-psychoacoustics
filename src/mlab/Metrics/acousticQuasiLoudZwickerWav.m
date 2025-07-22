function [loudness, fm, fn] = acousticQuasiLoudZwickerWav(p, sampleRatein, timeStep, axisN, fieldType, adjustEQL, ecmaEar)
% loudLevel = acousticQuasiLoudZwickerWav((p, sampleRatein, timeStep, axisN,
%                                         fieldType, adjustEQL, ecmaEar)
%
% Returns quasi-loudness using spectral elements of ISO 532-1 Zwicker
% loudness model for input pressure time-series. The temporal processing of
% the Zwicker model is disregarded. The input signal must be calibrated
% acoustic pressure (Pa).
%
% Optional modifications available comprise an adjustment to the spectral
% levels to improve agreement with the 2023 ISO 226 equal loudness
% contours, or (alternatively) the application of the ECMA-418-2:2024
% outer-middle ear filter responses (in 1/3-octaves) instead of the
% ISO 532-1:2017 outer ear transmission.
%
% Inputs
% ------
% p : vector or 2D matrix
%     the input signal as single mono or stereo audio (sound
%     pressure) signals
%
% sampleRatein : integer
%                the sample rate (frequency) of the input signal(s).
%
% timeStep : number
%            the time step value used to calculate the time-dependent Leq
%            and sharpness.
%
% axisN : integer (1 or 2, default: 1)
%         the time axis along which to calculate the sharpness.
%
% fieldType : keyword string (default: 'free-frontal')
%             determines whether the 'free-frontal' or 'diffuse' field stages
%             are applied in the outer-middle ear filter, or 'noOuter'
%             omits this filtering stage.
%
% adjustEQL : Boolean (default: false)
%             flag to indicate whether to apply adjustments for differences
%             between 1987 ISO 226 equal-loudness contours (which ISO 532-1
%             models) and 2023 ISO 226 equal-loudness contours. This option
%             is overridden if the ecmaEar option is set to true.
%
% ecmaEar : Boolean (default: false)
%           flag to indicate whether to substitute the outer-middle ear
%           filter response from ECMA-418-2:2024 for the ISO 532-1:2017
%           critical band ear transmission a0. This option overrides the
%           adjustEQL option.
%
% Returns
% -------
% fm : vector
%      the fractional octave band exact mid-frequencies for the input
%      spectra
%
% fn : vector
%      the fractional octave band nominal mid-frequencies for the input
%      spectra
%
% loudness : structure
%            contains the loudness output
%
% loudness contains the following outputs:
%
% loudTDep :  vector or matrix
%             time-dependent loudness
%             arranged as [time(, channels)]
% 
% loudPowAvg : number or vector
%              time-power-averaged loudness
%              arranged as [loudness(, channels)]
%
% loud5pcEx : number or vector
%             95th percentile (5% exceeded) loudness
%             arranged as [loudness(, channels)]
%
% loudLevel :  vector or matrix
%             time-dependent loudness level
%             arranged as [time(, channels)]
%
% loudLvlPowAvg : number or vector
%                 time-power-averaged loudness level
%                 arranged as [loudness(, channels)]
%
% specLoudAvg : vector or matrix
%               time-power-averaged specific loudness
%               arranged as [bands(, channels)]
%
% specLoudAvg : matrix
%               time-dependent specific loudness
%               arranged as [time, bands(, channels)]
%
% barkAxis : vector
%            critical band rates for specific loudness (0.1 dz intervals)
%
% timeOut : vector
%           time (seconds) corresponding with time-dependent outputs
%
% Assumptions
% -----------
% The input is a 1/3-octave band unweighted sound level spectrum or
% series of sound level spectra in dB re 2e-5 Pa.
%
% Requirements
% ------------
% None
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 23/04/2025
% Date last modified: 22/07/2025
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
        sampleRatein (1, 1) double {mustBePositive, mustBeInteger}
        timeStep (1, 1) double {mustBePositive}
        axisN (1, 1) {mustBeInteger, mustBeInRange(axisN, 1, 2)} = 1
        fieldType (1, :) string {mustBeMember(fieldType,...
                                                       {'free-frontal',...
                                                        'diffuse',...
                                                        'noOuter'})} = 'free-frontal'
        adjustEQL (1, 1) {mustBeNumericOrLogical} = false
        ecmaEar (1, 1) {mustBeNumericOrLogical} = false
    end

%% Input checks

% Orientate signal
if axisN == 2
    p = p.';
end

% Ensure sampling rate is full range
resampledRate = 48e3;
if sampleRatein ~= resampledRate  % Resample signal
    up = resampledRate/gcd(resampledRate, sampleRatein);  % upsampling factor
    down = sampleRatein/gcd(resampledRate, sampleRatein);  % downsampling factor
    p_re = resample(signal, up, down);  % apply resampling
else  % don't resample
    p_re = p;
end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Signal processing

% get number of channels
numChans = size(p_re, 2);

% Get time-averaged power spectrum
[pxx, ~] = poctave(p_re, resampledRate, 'spectrogram', 'BandsPerOctave', 3,...
                    'WindowLength', resampledRate*timeStep,...
                    'FrequencyLimits', [25, 12500]);

% reorientate power spectrum
if numChans == 1
    pxx = pxx.';
else
    pxx = permute(pxx, [2, 1, 3]);
end

% Calculate Leq
Leq = 10*log10(pxx/4e-10);

[loudness, fm, fn] = acousticLoudQuasiZwicker(Leq, [25, 12500], 2, fieldType, adjustEQL, ecmaEar);
loudness.timeOut = linspace(0, size(loudness.loudTDep, 1)*timeStep - timeStep, size(loudness.loudTDep, 1)).';

end  % end of acousticQuasiLoudZwickerWav function