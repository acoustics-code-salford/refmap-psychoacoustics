function [loudness, sharpness] = acousticSharpQuasiLoudWav(p, sampleRateIn, adjustLoud, timeStep, axisN, soundField, sharpMethod, outPlot, binaural, recalibrate)
% [loudness, sharpness] = acousticSharpQuasiLoudWav(p, sampleRatein, adjustLoud,
%                                                   timeStep, axisN,
%                                                   soundField, sharpMethod,
%                                                   outPlot, binaural, recalibrate)
%
% Returns quasi-sharpness and quasi-sharpness values using
% acousticSharpFromQuasi.m from quasi-loudness results obtained using
% acousticQuasiLoudZwicker.m with Leq spectra obtained with poctave.
%
% The sharpness model used can be specified using the 'sharpMethod' input
% argument. Options comprise 'aures', 'vonbismarck', or 'widmann' (which
% is the model standardised in DIN 45692:2009).
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
% adjustLoud : keyword string
%   Indicates whether to apply adjustments for the outer-middle ear filter
%   response and loudness transformations from ISO 532-3:2023 ('iso5323')
%   or ECMA-418-2:2024 ('ecma4182').
%   The default option ('iso5321') applies no adjustment, so follows ISO
%   532-1 more closely (for closer agreement with Zwicker's model).
%
% timeStep : number
%   Time step value (seconds) used to calculate the time-dependent Leq
%   and sharpness.
%
% axisN : integer (1 or 2, default: 1)
%   Time axis along which to calculate the sharpness.
%
% soundField : keyword string (default: 'freeFrontal')
%   Determines whether the 'freeFrontal' or 'diffuse' field stages
%   are applied in the outer/outer-middle ear filter. The
%   'noOuter' option also allows either no transmission filter
%   (if adjustLoud is not 'ecma4182'), or omits the outer ear filter stage
%   from the ECMA-418-2:2024 filter (if adjustLoud is 'ecma4182').
%
% sharpMethod : keyword string (default: 'aures')
%   The sharpness weighting method to apply. Options: 'aures',
%   'vonbismarck', 'widmann'.
%
% outPlot : Boolean true/false (default: false)
%   Flag indicating whether to generate a figure from the output.
% 
% binaural : Boolean true/false (default: true)
%   flag indicating whether to output binaural sharpness for
%   stereo input signal. (It is assumed the relationship for binaural
%   sharpness follows that of binaural loudness, which seems to be
%   supported by available evidence https://doi.org/10.1051/aacus/2025048)
%
% recalibrate : Boolean (default: false)
%   Indicates whether to apply a model-specific calibration to predict
%   absolute loudness values (derived from empirical data).
%
% Returns
% -------
%
% sharpness : structure
%   contains the output.
%
% sharpness contains the following outputs:
%
% sharpnessTDep : vector or matrix
%   Time-dependent sharpness arranged as [time(, channels)].
% 
% sharpnessPowAvg : number or vector
%   Time-power-averaged sharpness arranged as [sharpness(, channels)].
%
% sharpness5pcEx : number or vector
%   95th percentile (5% exceeded) sharpness
%   arranged as [sharpness(, channels)].
%
% soundField : string
%   Identifies the soundfield type applied (the input argument soundField)
%
% timeOut : vector
%   Time (seconds) corresponding with time-dependent outputs.
%
% sharpWeight : string
%   Indicates which sharpness weighting method was applied.
%
% adjustLoud : string
%   Indicates which loudness model standard was adjusted for (=adjustLoud
%   input).
%
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
% If outPlot=true, plots are returned illustrating the time-dependent
% loudness and sharpness, along with the time-aggregated values.
% A plot is returned for each input channel.
%
% Assumptions
% -----------
% The input signal is calibrated to units of acoustic pressure in Pascals
% (Pa).
%
% References
% ----------
% Aures model is described in:
%
% Aures, W., 1985. Berechnungsverfahren für den sensorischen Wohlklang
% beliebiger Schallsignale, Acta Acustica united with Acustica, 59(2),
% 130-141.
% https://www.ingentaconnect.com/content/dav/aaua/1985/00000059/00000002/art00008
%
% von Bismarck-Zwicker model is described in:
% 
% von Bismarck, G., 1974. Sharpness as an attribute of the timbre of
% steady sounds, Acta Acustica united with Acustica, 30(3) 159-172.
% https://www.ingentaconnect.com/content/dav/aaua/1974/00000030/00000003/art00006
%
% and:
%
% Zwicker, E. & Fastl, H., 1999. Psychoacoustics: Facts and models.
% Springer-Verlag.
%
% Widmann model is described in DIN 45692:2009
%
% Zwicker loudness is defined in ISO 532-1:2017. Modifications are based on
% ISO 532-3:2023 and ECMA-418-2:2025.
%
% The assumed binaural perception of sharpness is based on evidence found
% in:
%
% Hochbaum, F, Hundt,  T, Fiebig, A & Brinkmann, F, 2025. Directional
% sharpness perception under different listening conditions, Acta Acustica,
% 9, 60. https://doi.org/10.1051/aacus/2025048
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
% Date created: 30/04/2025
% Date last modified: 08/03/2026
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
% This code calls sub-component file 'cmap_viridis.txt'. The contents of
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
        adjustLoud (1, :) string {mustBeMember(adjustLoud,...
                                               {'iso5321',...
                                                'iso5323',...
                                                'ecma4182'})} = 'iso5321'
        timeStep (1, 1) double {mustBePositive} = 0.1
        axisN (1, 1) {mustBeInteger, mustBeInRange(axisN, 1, 2)} = 1
        soundField (1, :) string {mustBeMember(soundField,...
                                              {'freeFrontal',...
                                               'diffuse',...
                                               'noOuter'})} = 'freeFrontal'
        sharpMethod (1, :) string {mustBeMember(sharpMethod,...
                                                {'aures',...
                                                 'vonbismarck',...
                                                 'widmann'})} = 'aures'
        outPlot (1, 1) {mustBeNumericOrLogical} = false
        binaural (1, 1) {mustBeNumericOrLogical} = true
        recalibrate (1, 1) {mustBeNumericOrLogical} = true
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Signal processing

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

% get number of channels
numChans = size(p_re, 2);

% Get time-averaged power spectrum
[pxx, ~] = poctave(p_re, sampleRateIn, 'spectrogram', 'BandsPerOctave', 3,...
                   'WindowLength', sampleRateIn*timeStep,...
                   'FrequencyLimits', [25, 12500]);

% reorientate power spectrum
if numChans >= 2
    pxx = permute(pxx, [2, 1, 3]);
else
    pxx = pxx.';
end

% Calculate Leq
Leq = 10*log10(pxx/4e-10);

% Calculate loudness
loudness = acousticQuasiLoudZwicker(Leq, [25 12600], timeStep, 2, soundField,...
                                    adjustLoud, outPlot, recalibrate);

% Calculate sharpness
sharpness = acousticSharpFromQuasiLoud(loudness.loudTDep,...
                                       loudness.specLoud, adjustLoud,...
                                       [25 12600], timeStep, sharpMethod,...
                                       outPlot, binaural);

end  % end of acousticSharpQuasiLoudWav function
