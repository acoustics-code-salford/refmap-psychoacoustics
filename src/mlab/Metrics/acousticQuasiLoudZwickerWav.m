function loudness = acousticQuasiLoudZwickerWav(p, sampleRateIn, timeStep, axisN, soundField, adjustLoud)
% loudness = acousticQuasiLoudZwickerWav(p, sampleRateIn, timeStep, axisN,
%                                        soundField, adjustLoud)
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
%   The input signal as single mono or stereo audio (sound
%   pressure) signals.
%
% sampleRateIn : integer
%   The sample rate (frequency) of the input signal(s).
%
% timeStep : number
%   The time step value used to calculate the time-dependent Leq
%   and loudness.
%
% axisN : integer (1 or 2, default: 1)
%   The time axis along which to calculate the loudness.
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
% Returns
% -------
% loudness : structure
%   contains the loudness output
%
% loudness contains the following outputs:
%
% loudTDep :  vector or matrix
%   time-dependent loudness
%   arranged as [time(, channels)]
% 
% loudPowAvg : number or vector
%   time-power-averaged loudness
%   arranged as [loudness(, channels)]
%
% loud5pcEx : number or vector
%   95th percentile (5% exceeded) loudness
%   arranged as [loudness(, channels)]
%
% loudLevel :  vector or matrix
%   time-dependent loudness level
%   arranged as [time(, channels)]
%
% loudLvlPowAvg : number or vector
%   time-power-averaged loudness level
%   arranged as [loudness(, channels)]
%
% specLoudAvg : vector or matrix
%   time-power-averaged specific loudness
%   arranged as [bands(, channels)]
%
% specLoudAvg : matrix
%   time-dependent specific loudness
%   arranged as [time, bands(, channels)]
%
% barkAxis : vector
%   critical band rates for specific loudness (0.1 dz intervals)
%
% timeOut : vector
%   time (seconds) corresponding with time-dependent outputs
%
% freqInMid : vector
%   The 1/3-octave band exact mid-frequencies for the input spectra.
%
% freqInNom : vector
%   The 1/3-octave band nominal mid-frequencies for the input spectra.
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
% Date last modified: 05/12/2025
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

% get number of channels
chansIn = size(p_re, 2);

% Get time-averaged power spectrum
[pxx, ~] = poctave(p_re, resampledRate, 'spectrogram', 'BandsPerOctave', 3,...
                   'FilterOrder', 6,...
                   'FrequencyLimits', [25, 12600],...
                   'Weighting', 'none',...
                   'WindowLength', resampledRate*timeStep,...
                   'OverlapPercent', 0);

% reorientate power spectrum
if chansIn == 1
    pxx = pxx.';
else
    pxx = permute(pxx, [2, 1, 3]);
end

% Calculate Leq
Leq = 10*log10(pxx/4e-10);

loudness = acousticQuasiLoudZwicker(Leq, [25, 12500], 2,...
                                    soundField, adjustLoud);
loudness.timeOut = linspace(0, size(loudness.loudTDep, 1)*timeStep...
                            - timeStep, size(loudness.loudTDep, 1)).';

end  % end of acousticQuasiLoudZwickerWav function