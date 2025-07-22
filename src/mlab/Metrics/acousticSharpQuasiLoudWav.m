function sharpnessQZ = acousticSharpQuasiLoudWav(p, sampleRatein, timeStep, axisN, fieldType, sharpMethod, adjustEQL, ecmaEar, outPlot)
% sharpnessQZ = acousticSharpQuasiLoudWav(p, sampleRatein, timeStep, axisN,
%                                     fieldType, sharpMethod,
%                                     adjustEQL, ecmaEar, outPlot)
%
% Returns quasi-sharpness values using acousticSharpFromQuasiZwicker.m from 
% quasi-loudness results obtained using acousticQuasiLoudZwicker.m with Leq spectra
% obtained with poctave.
%
% The sharpness model used can be specified using the 'method' input
% argument. Options comprise 'aures', 'vonBismarck', or 'widmann' (which
% is the model standardised in DIN 45692:2009).
%
% Optional modifications to the loudness calculation comprise an adjustment
% to the spectral levels to improve agreement with the 2023 ISO 226 equal
% loudness contours, or (alternatively) the application of the ECMA-418-2:2024
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
%             are applied in the outer/outer-middle ear filter. The
%             'noOuter' option also allows either no transmission filter
%             (if ecmaEar is false), or omits the outer ear filter stage
%             from the ECMA-418-2:2024 filter (if ecmaEar is true).
%
% sharpMethod : keyword string (default: 'aures')
%               the sharpness method to apply. Options: 'aures', 'vonbismarck',
%               'widmann'.
%
% outplot : Boolean true/false (default: false)
%           flag indicating whether to generate a figure from the output
%
% binaural : Boolean true/false (default: true)
%            flag indicating whether to output binaural sharpness for
%            stereo input signal. (Experimental: it is assumed the
%            relationship for binaural sharpness follows that of binaural
%            SHM loudness.
% 
% Returns
% -------
%
% sharpnessQZ : structure
%                     contains the output
%
% sharpnessQZ contains the following outputs:
%
% sharpnessTDep : vector or matrix
%                 time-dependent sharpness
%                 arranged as [time(, channels)]
% 
% sharpnessPowAvg : number or vector
%                   time-power-averaged sharpness
%                   arranged as [sharpness(, channels)]
%
% sharpness5pcEx : number or vector
%                  95th percentile (5% exceeded) sharpness
%                  arranged as [sharpness(, channels)]
%
% soundField : string
%              identifies the soundfield type applied (the input argument
%              fieldtype)
%
% timeOut : vector
%           time (seconds) corresponding with time-dependent outputs
%
% method : string
%          indicates which sharpness method was applied
%
% If outplot=true, a plot is returned illustrating the
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
% Requirements
% ------------
% None
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 30/04/2025
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
        fieldType (1, :) string {mustBeMember(fieldType,...
                                              {'free-frontal',...
                                               'diffuse',...
                                               'noOuter'})} = 'free-frontal'
        sharpMethod (1, :) string {mustBeMember(sharpMethod,...
                                                {'aures',...
                                                 'vonbismarck',...
                                                 'widmann'})} = 'aures'
        adjustEQL (1, 1) {mustBeNumericOrLogical} = false
        ecmaEar (1, 1) {mustBeNumericOrLogical} = false
        outPlot {mustBeNumericOrLogical} = false
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
[loudness, ~, ~] = acousticLoudQuasiZwicker(Leq, [25 12500], 2, fieldType, adjustEQL, ecmaEar);

% Calculate sharpness
sharpnessQZ = acousticSharpFromQuasiZwicker(loudness.loudTDep, loudness.specLoud, timeStep, sharpMethod, outPlot);

end  % end of acousticSharpQuasiLoudWav function
