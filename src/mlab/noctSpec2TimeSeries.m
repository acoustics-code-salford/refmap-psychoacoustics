function signal = noctSpec2TimeSeries(specRMSdB, dBRef, fLim, bandsPerOctave, filterOrder, duration, sampleRate)
% signal = noctSpec2TimeSeries(specdB, dBRef, fLim, bandsPerOctave, filterOrder, sRateOut)
%
% Returns
%
% Inputs
% ------
% specdB : vector
%   Input RMS spectrum in decibels referenced to input argument dBRef.
%
% Returns
% -------
% signal : vector
%   Output signal in units corresponding with the dB input reference.
%
% Assumptions
% -----------
%
% Requirements
% ------------
% Signal Processing Toolbox
% Audio Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Authors: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 03/12/2025
% Date last modified: 03/12/2025
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
%% Arguments validation

% fractional octave frequencies
[freqMid, ~, ~, ~] = noctf(fLim, bandsPerOctave);

%% Generate initial white noise
N  = round(sampleRate*duration);
x  = randn(N, 1);

switch bandsPerOctave
    case 1
        bandwidth = '1 octave';
    case 2
        bandwidth = '1/2 octave';
    case 3
        bandwidth = '1/3 octave';
    case 6
        bandwidth = '1/6 octave';
    case 12
        bandwidth = '1/12 octave';
    case 24
        bandwidth = '1/24 octave';
    case 48
        bandwidth = '1/48 octave';
end


%% Design filters
for kk = numel(freqMid):-1:1
    filters{kk} = octaveFilter('FilterOrder', filterOrder, ...
                               'CenterFrequency', freqMid(kk), ...
                               'Bandwidth', bandwidth, ...
                               'SampleRate', sampleRate);
end

%% Filter, scale each band, and sum
signal = zeros(N,1);
for kk = 1:numel(freqMid)
    % Band-passed noise
    bandNoise = filters{kk}(x);

    % Current RMS
    bandNoiseRMS = rms(bandNoise);

    % Target RMS
    bandNoiseTargetRMS = db2mag(specRMSdB(kk))*dBRef;

    % Scale and sum
    signal = signal + bandNoise*(bandNoiseTargetRMS/bandNoiseRMS);
end
