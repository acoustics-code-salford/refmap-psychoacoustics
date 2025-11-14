function sharpness = acousticSharpQuasiLoudWav(p, sampleRatein, timeStep, axisN, soundField, sharpMethod, adjustLoud, outPlot, binaural)
% sharpness = acousticSharpQuasiLoudWav(p, sampleRatein, timeStep, axisN,
%                                       soundField, sharpMethod,
%                                       adjustLoud, outPlot, binaural)
%
% Returns quasi-sharpness values using acousticSharpFromQuasi.m from 
% quasi-loudness results obtained using acousticQuasiLoudZwicker.m with Leq spectra
% obtained with poctave.
%
% The sharpness model used can be specified using the 'sharpMethod' input
% argument. Options comprise 'aures', 'vonbismarck', or 'widmann' (which
% is the model standardised in DIN 45692:2009).
%
% Optional modifications available comprise an adjustment to the spectral
% levels to improve agreement with the 2023 ISO 226 equal loudness
% contours ('iso226'), or (alternatively) the application of the
% ECMA-418-2:2025 outer-middle ear filter responses (in 1/3-octaves)
% instead of the ISO 532-1:2017 outer ear transmission, and an
% approximated non-linear transformation following ECMA-418-2:2025
% ('ecma4182'), instead of Zwicker's original.
%
% Inputs
% ------
% p : vector or 2D matrix
%   the input signal as single mono or stereo audio (sound
%   pressure) signals
%
% sampleRatein : integer
%   the sample rate (frequency) of the input signal(s).
%
% timeStep : number
%   the time step value used to calculate the time-dependent Leq
%   and sharpness.
%
% axisN : integer (1 or 2, default: 1)
%   the time axis along which to calculate the sharpness.
%
% soundField : keyword string (default: 'freeFrontal')
%   determines whether the 'freeFrontal' or 'diffuse' field stages
%   are applied in the outer/outer-middle ear filter. The
%   'noOuter' option also allows either no transmission filter
%   (if adjustLoud is not 'ecma4182'), or omits the outer ear filter stage
%   from the ECMA-418-2:2024 filter (if adjustLoud is 'ecma4182').
%
% sharpMethod : keyword string (default: 'aures')
%   the sharpness method to apply. Options: 'aures', 'vonbismarck',
%   'widmann'.
%
% adjustLoud : keyword string (default: 'none')
%   indicates whether to apply adjustments for:
%   ('iso226') the differences between 1987 ISO 226 equal-loudness contours
%   (which ISO 532-1 models) and 2023 ISO 226 equal-loudness contours, or;
%   ('ecma4182') the outer-middle ear filter response from ECMA-418-2:2024
%   omitting the ISO 532-1:2017 critical band ear transmission a0, and
%   adapting the loudness transformation to agree more closely with the
%   ECMA-418-2:2024 transformation.
%   The default option ('none') applies no adjustment, so follows ISO 532-1
%   more closely (for closer agreement with Zwicker's model).
%
% outplot : Boolean true/false (default: false)
%   flag indicating whether to generate a figure from the output
% 
% binaural : Boolean true/false (default: true)
%   flag indicating whether to output binaural sharpness for
%   stereo input signal. (It is assumed the relationship for binaural
%   sharpness follows that of binaural loudness, which seems to be
%   supported by available evidence https://doi.org/10.1051/aacus/2025048)
%
% Returns
% -------
%
% sharpness : structure
%   contains the output
%
% sharpness contains the following outputs:
%
% sharpnessTDep : vector or matrix
%   time-dependent sharpness arranged as [time(, channels)]
% 
% sharpnessPowAvg : number or vector
%   time-power-averaged sharpness arranged as [sharpness(, channels)]
%
% sharpness5pcEx : number or vector
%   95th percentile (5% exceeded) sharpness
%   arranged as [sharpness(, channels)]
%
% soundField : string
%   identifies the soundfield type applied (the input argument fieldtype)
%
% timeOut : vector
%   time (seconds) corresponding with time-dependent outputs
%
% method : string
%   indicates which sharpness method was applied
%
% If outPlot=true, a plot is returned illustrating the
% time-dependent sharpness, with the time-aggregated values.
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
% von Bismarck model is described in:
% 
% von Bismarck, G., 1974. Sharpness as an attribute of the timbre of
% steady sounds, Acta Acustica united with Acustica, 30(3) 159-172.
% https://www.ingentaconnect.com/content/dav/aaua/1974/00000030/00000003/art00006
%
% Widmann model is described in DIN 45692:2009
%
% Zwicker loudness is defined in ISO 532-1:2017. Modifications are based on
% ISO 226:2023, ISO 226:1987 and ECMA-418-2:2024.
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
% Date last modified: 13/11/2025
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
        sampleRatein (1, 1) double {mustBePositive, mustBeInteger}
        timeStep (1, 1) double {mustBePositive}
        axisN (1, 1) {mustBeInteger, mustBeInRange(axisN, 1, 2)} = 1
        soundField (1, :) string {mustBeMember(soundField,...
                                              {'freeFrontal',...
                                               'diffuse',...
                                               'noOuter'})} = 'freeFrontal'
        sharpMethod (1, :) string {mustBeMember(sharpMethod,...
                                                {'aures',...
                                                 'vonbismarck',...
                                                 'widmann'})} = 'aures'
        adjustLoud (1, :) string {mustBeMember(adjustLoud,...
                                               {'none',...
                                                'iso226',...
                                                'ecma4182'})} = 'none'
        outPlot {mustBeNumericOrLogical} = false
        binaural {mustBeNumericOrLogical} = true
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Signal processing

% Orientate signal
if axisN == 2
    p = p.';
end

% get number of channels
numChans = size(p, 2);

% Get time-averaged power spectrum
[pxx, ~] = poctave(p, sampleRatein, 'spectrogram', 'BandsPerOctave', 3,...
                   'WindowLength', sampleRatein*timeStep,...
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
loudness = acousticQuasiLoudZwicker(Leq, [25 12600], 2, soundField,...
                                    adjustLoud);

% Calculate sharpness
sharpness = acousticSharpFromQuasiLoud(loudness.loudTDep,...
                                       loudness.specLoud, adjustLoud,...
                                       timeStep, sharpMethod,...
                                       outPlot, binaural);

end  % end of acousticSharpQuasiLoudWav function
