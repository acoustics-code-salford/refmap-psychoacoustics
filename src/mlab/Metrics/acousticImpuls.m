function impulsive = acousticImpuls(p, sampleRateIn, axisN, boundSkip, soundField, loudMethod, outPlot)
% impulsive = acousticImpuls(p, sampleRateIn, axisN, startSkip,
%                            soundField, loudMethod, outPlot)
%
% Returns impulsive loudness values according to Willemsen & Rao (2010) for
% an input calibrated single mono or single stereo audio (sound pressure)
% time-series signal, p.
%
% p : vector or 2D matrix
%   Input signal as single mono or stereo audio (sound pressure) signals.
%
% sampleRateIn : integer
%   Sample rate (frequency) of the input signal(s).
%
% axisN : integer (1 or 2, default: 1)
%   Time axis along which to calculate the tonality.
%
% boundSkip : number (default: 0.5)
%   Amount of time to skip at the start of the signal for
%   calculating time-aggregated outputs (starts from next input
%   sample), to account for uncertain boundary behaviour (loudness
%   calculation filters and percentile windowing). The default is the
%   minimum value allowed and is skipped from both start and end
%   boundaries.
%
% soundField : keyword string (default: 'freeFrontal')
%   Determines whether the 'freeFrontal' or 'diffuse' field stages
%   are applied in the outer-middle ear filtering.
%
% loudMethod : keyword string (default: 'zwicker')
%   Determines the loudness model to use (options: 'zwicker' or 'sottek').
%   If 'zwicker', ISO 532-1:2017 loudness is used; if 'sottek',
%   ECMA-418-2:2025 loudness is used. The Sottek Hearing Model loudness is
%   a more recent model that agrees better with modern equal loudness
%   contours, enables combined binaural outputs (for 2-channel inputs),
%   and is supported by listening test data for a wide range of
%   sounds. The adapted Zwicker loudness model is older, but is less
%   computationally intensive, and features a higher output sampling rate,
%   which may be more suited for characterising highly impulsive sounds
%   with very rapid transients. Note that if 'sottek' is used with a
%   2-channel input, the single-channel output corresponds with the
%   combined binaural loudness.
%
% outPlot : Boolean true/false (default: false)
%   Flag indicating whether to generate a figure from the output.
%
% Returns
% -------
% impulsive : structure
%   contains the output
%
% impulsive contains the following outputs:
%
% impulsLoudTDep : vector or matrix
%   Time-dependent impulsive loudness arranged as [time(, channels)].
% 
% impulsLoudAvg : number or vector
%   Time-averaged impulsive loudness arranged as [sharpness(, channels)].
%
% impulsLoudPowAvg : number or vector
%   Time-power-averaged impulsive loudness arranged as [sharpness(, channels)].
%
% timeOut : vector
%   Time (seconds) corresponding with time-dependent outputs.
%
% soundField : string
%              identifies the soundfield type applied (= input argument)
%
% If outplot=true, a set of plots is returned illustrating the
% time-dependent impulsiveness, with the time-aggregated values.
% A set of plots is returned for each input channel.
%
% Assumptions
% -----------
% The input signal is calibrated to units of acoustic pressure in Pascals
% (Pa).
%
% References
% ----------
% The impulsiveness model is described in:
%
% Willemsen, A M & Rao, M D, 2010. Characterization of sound quality of
% impulsive sounds using loudness based metric. In: Proceedings of 20th
% International Congress on Acoustics (ICA) 2010, 23–27 August 2010,
% Sydney, Australia.
%
% Requirements
% ------------
% Audio Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 18/11/2025
% Date last modified: 17/12/2025
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
        axisN (1, 1) {mustBeInteger, mustBeInRange(axisN, 1, 2)} = 1
        boundSkip (1, 1) {mustBeGreaterThanOrEqual(boundSkip, 0.5)} = 0.5    
        soundField (1, :) string {mustBeMember(soundField,...
                                               {'freeFrontal',...
                                                'diffuse'})} = 'freeFrontal'
        loudMethod (1, :) string {mustBeMember(loudMethod,...
                                               {'zwicker',...
                                                'sottek'})} = 'zwicker'
        outPlot (1, 1) {mustBeNumericOrLogical} = false
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Input checks

% Orient input matrix
if axisN == 2
    p = p.';
end

% Check the length of the input data (must be longer than 300 ms)
if size(p, 1) <=  (2*boundSkip)*sampleRateIn
    error('Error: Input signal is too short along the specified axis to calculate sharpness (must be longer than 2*boundSkip s)')
end

% Check the channel number of the input data
if size(p, 2) > 2
    error('Error: Input signal comprises more than two channels')
else
    chansIn = size(p, 2);
    if chansIn > 1
        chans = ["Stereo left";
                 "Stereo right"];
    else
        chans = "Mono";
    end
end

%% Signal processing

if strcmp(loudMethod, 'zwicker')
    % calibration value to ensure reference signal equals 1 sone (reference:
    % sinusoid at 1 kHz at 40 dB free-field)
    calN = 0.999511042933208;
    % convert soundField keyword for compatibility with Audio Toolbox
    if strcmp(soundField, 'freeFrontal')
        soundFieldKey = "free";
    else
        soundFieldKey = "diffuse";
    end
    
    % calculate loudness
    [loudness, ~, ~] = acousticLoudness(p, sampleRateIn, 1,...
                                        'Method', 'ISO 532-1',...
                                        'TimeVarying', true,...
                                        'SoundField', soundFieldKey,...
                                        'TimeResolution', 'high');
    dt = 5e-4;
    sampleRateOut = 1/dt;

    % calibrate time-dependent loudness to ensure 1 sone for ref signal
    loudnessTDep = calN*loudness;

else
    % calculate loudness
    loudnessSHM = acousticSHMLoudness(p, sampleRateIn, axisN,...
                                      soundField, false, false, true);

    sampleRateOut = 187.5;
    dt = 1/sampleRateOut;

    % extract binaural or monaural loudness
    if chansIn == 2
        loudnessTDep = loudnessSHM.loudnessTDepBin;
    else
        loudnessTDep = loudnessSHM.loudnessTDep;
    end
end

% output time axis
timeOut = transpose(0:dt:size(loudnessTDep, 1)/sampleRateOut - dt);
sampleN = numel(timeOut);

% calculate background/baseline loudness
% indices for calculation
startIdx = find(timeOut >= 0.5, 1, 'first');
endIdx = find(timeOut <= timeOut(end) - 0.5, 1, 'last');
loudnessBase = loudnessTDep;

% start section
for it = 1:startIdx - 1
    startWin = find(timeOut >= 0, 1, 'first');
    endWin = find(timeOut <= timeOut(it) + 0.5, 1, 'last');
    loudnessBase(it,:) = prctile(loudnessTDep(startWin:endWin ,:), 5, 1);
end

% end section
for it = endIdx + 1:sampleN
    startWin = find(timeOut >= timeOut(it) - 0.5, 1, 'first');
    endWin = sampleN;
    loudnessBase(it,:) = prctile(loudnessTDep(startWin:endWin ,:), 5, 1);
end

% remaining samples
startWin = find(timeOut >= (timeOut(startIdx) - 0.5), 1, 'first');
endWin = find(timeOut <= (timeOut(startIdx) + 0.5), 1, 'last');
for it = startIdx:endIdx
    loudnessBase(it, :) = prctile(loudnessTDep(startWin:endWin, :), 5, 1);
    startWin = startWin + 1;
    endWin = endWin + 1;
end

% impulsive loudness (loudness difference)
impulsLoudTDep = max(loudnessTDep - loudnessBase, 0);

% aggregation indices
aggStartIdx = find(timeOut > boundSkip, 1, 'first');
aggEndIdx = find(timeOut < timeOut(end) - boundSkip, 1, 'last');

% impulsive loudness average
impulsLoudAvg = mean(impulsLoudTDep(aggStartIdx:aggEndIdx, :), 1);

% impulsive loudness power average
impulsLoudPowAvg = (sum(impulsLoudTDep(aggStartIdx:aggEndIdx, :).^(1/log10(2)), 1)./size(impulsLoudTDep(aggStartIdx:aggEndIdx, :), 1)).^(log10(2));

%% Output plotting
% TODO
if outPlot
    disp(chans)
end

%% Assign outputs

impulsive.impulsLoudTDep = impulsLoudTDep;
impulsive.impulsLoudAvg = impulsLoudAvg;
impulsive.impulsLoudPowAvg = impulsLoudPowAvg;
impulsive.timeOut = timeOut;
impulsive.soundField = soundField;
