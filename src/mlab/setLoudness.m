function [pSet, soneOut] = setLoudness(p, sampleRateIn, soneTarget, axisN, soundField, aggregate, method, binaural, precision, timeLimitX)
% pSet = setLoudness(p, sampleRateIn, soneTarget, axisN, soundField, aggregate, method, precision, timeLimitX)
%
% Returns an output sound pressure time-series signal 'pSet' calibrated to
% the target loudness, 'soneTarget', for an input calibrated single mono or
% single stereo audio sound pressure time-series signal, 'p'. The loudness model
% is set with the input 'method'. The value of soneTarget corresponds with the time
% aggregation approach defined using 'aggregate'. For a stereo input, if
% 'binaural' is true and the method is ECMA-418-2, the quadratic mean
% estimation of combined binaural loudness will be applied.
% Otherwise, the loudness value will be calibrated for the channel with the
% highest loudness.
%
% Inputs
% ------
% p : vector or 2D matrix
%     the input signal as single mono or stereo audio (sound
%     pressure) signals
%
% sampleRateIn : integer
%                the sample rate (frequency) of the input signal(s)
%
% soneIn : float
%          the desired (target) loudness, in sones
%
% axisN : integer (1 or 2, default: 1)
%         the time axis along which to calculate the loudness
%
% soundField : keyword string (default: 'freeFrontal')
%              determines whether the 'freeFrontal' or 'diffuse' field stages
%              are applied in the outer-middle ear filtering.
%
% aggregate : keyword string (default: PowAvg)
%             the time aggregation corresponding with the value of sone.
%             Options include 'PowAvg' (ECMA-418-2:2025), 'Mean', 'Max',
%             'Min', or '5PcEx' (5%-exceeded).
%
% method : keyword string (default: ECMA-418-2)
%          the loudness model to use - options include ECMA-418-2:2025
%          ('ECMA-418-2', default) or ISO 532-1:2017 ('ISO-532-1').
%
% binaural : Boolean (default: true)
%            flag to indicate whether to take the combined binaural
%            loudness as the target in the case of a stereo input (true),
%            or otherwise use the loudness for the channel with the higher
%            loudness.
%
% precision : float (default: 0.05)
%             the target precision within which the algorithm will seek to
%             match the target loudness (> 0, < 1).
%
% timeLimitX : float (default: 10)
%              a multiplier that defines a time limit for the algorithm to
%              attempt loudness matching within, based on the input signal.
%              The default time is 10x the input signal duration (i.e., for
%              an input signal of 10 s duration, the algorithm will stop
%              trying to match loudness after 100 s, and will return output
%              based on the last nearest value).
%
% Returns
% -------
%
% pSet : vector or 2D matrix
%        the output calibrated signal as single mono or stereo audio (sound
%        pressure) signals
%
% soneOut : float
%          the achieved output loudness, in sones
%
% Assumptions
% -----------
% The input signal is calibrated to units of acoustic pressure in Pascals
% (Pa).
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
% Date created: 10/07/2025
% Date last modified: 11/07/2025
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
% This code calls sub-component file 'cmap_inferno.txt'. The contents of
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
        soneTarget (1, 1) {mustBePositive}
        axisN (1, 1) {mustBeInteger, mustBeInRange(axisN, 1, 2)} = 1
        soundField (1, :) string {mustBeMember(soundField,...
                                               {'freeFrontal',...
                                                'diffuse'})} = 'freeFrontal'
        aggregate (1, :) string {mustBeMember(aggregate,...
                                              {'PowAvg',...
                                               'Mean', 'Max', 'Min',...
                                               '5PcEx'})} = 'PowAvg'
        method (1, :) string {mustBeMember(method,...
                                           {'ECMA-418-2',...
                                            'ISO 532-1'})} = 'ECMA-418-2'
        binaural (1, 1) {mustBeNumericOrLogical} = true
        precision (1, 1) double {mustBePositive, mustBeLessThan(precision, 1)} = 0.05
        timeLimitX (1, 1) double {mustBePositive} = 10
    end

%% Load path
addpath(genpath(fullfile("refmap-psychoacoustics", "src", "mlab")))

%% Input checks
% Orient input matrix
if axisN == 2
    p = p.';
end

if strcmp(method, 'ECMA-418-2')
    startTSkip = 0.3;
else 
    startTSkip = 0.5;
end

% Check the length of the input data (must be longer than the min duration s)
if size(p, 1)/sampleRateIn <=  startTSkip
    error("Error: Input signal is too short along the specified axis to calculate loudness (must be longer than " + num2str(startTSkip, 1) + " s)")
end

% check channel number
if size(p, 2) > 2
    error('Error: Input signal comprises more than two channels')
else
    chansIn = size(p, 2);
end

%% Define constants

% convert soundField keyword for compatibility with Audio Toolbox
if strcmp(soundField, 'freeFrontal')
    soundFieldKey = "free";
else
    soundFieldKey = "diffuse";
end

% assign aggregation sample skip based on output sample rate
switch method
    case 'ECMA-418-2'
        startSkipIdx = ceil(startTSkip*187.5) + 1;
    otherwise
        startSkipIdx = ceil(startTSkip*500) + 1;
end

% determine timeLimit in seconds
timeLimit = size(p, 1)/sampleRateIn*timeLimitX;


%% Signal processing

if (~binaural && chansIn == 2) || strcmp(method, 'ISO 532-1')
    N = runLoudness(p, aggregate, method);
    [~, whichChan] = max(N, [], 2);
    % use only the higher loudness signal to calculate the loudness
    pTarget = p(:, whichChan);
else
    whichChan = 1;
    pTarget = p;
end

% get the loudness for the input signal
N = runLoudness(pTarget, aggregate, method);

% calculate the difference compared with the target value
loudDiff = abs(N - soneTarget);

% if loudness does not already meet target, loop over loudness calculation
% and iterate signal calibration to reduce loudness difference between
% signal and target
if loudDiff > precision
    
    tic  % start timer    
    iterIdx = 1;
    while toc <= timeLimit && loudDiff > precision      

        divVal = N/soneTarget;
        pTarget = pTarget/divVal;

        N = runLoudness(pTarget, aggregate, method);

        loudDiff = abs(N - soneTarget);
    
        disp("Iteration " + num2str(iterIdx) + " loudness: " + num2str(N) + " sone")
        disp("Iteration " + num2str(iterIdx) + " loudness difference from target: " + num2str(loudDiff) + " sone")

        iterIdx = iterIdx + 1;
    end

    if toc > timeLimit
        disp("Time limit of " + num2str(timeLimit) + " seconds exceeded! Outputting nearest value obtainted. Increase time limit to achieve precision.")
    elseif loudDiff <= precision
        disp("Target loudness achieved to input precision.")
    end

    pSet = p/(rms(p(:, whichChan))/rms(pTarget));
    if ~binaural && chansIn == 2
        soneOut = runLoudness(pSet, aggregate, method);
    else
        soneOut = N;
    end

else  % if precision already met
    pSet = p;
    if ~binaural && chansIn == 2
        soneOut = runLoudness(pSet, aggregate, method);
    else
        soneOut = N;
    end
end

    % nested function to run loudness calculations
    function N = runLoudness(signal, aggregate, method)
        switch method
            case 'ECMA-418-2'
                loudness = acousticSHMLoudness(signal, sampleRateIn, axisN, soundField, false, false, binaural);

                if binaural && chansIn == 2
                    if strcmp(aggregate, 'PowAvg')
                        N = loudness.loudnessPowAvgBin;
                    else
                        switch aggregate
                            case 'Mean'
                                N = mean(loudness.loudnessTDepBin(startSkipIdx:end));
                            case 'Max'
                                N = max(loudness.loudnessTDepBin(startSkipIdx:end));
                            case 'Min'
                                N = min(loudness.loudnessTDepBin(startSkipIdx:end));
                            case '5PcEx'
                                N = prctile(loudness.loudnessTDepBin(startSkipIdx:end), 95);
                        end  % end of switch for aggregate
                    end  % end of if branch for aggregate
                else
                    if strcmp(aggregate, 'PowAvg')
                        N = loudness.loudnessPowAvg;
                    else
                        switch aggregate
                            case 'Mean'
                                N = mean(loudness.loudnessTDep(startSkipIdx:end, :));
                            case 'Max'
                                N = max(loudness.loudnessTDep(startSkipIdx:end, :));
                            case 'Min'
                                N = min(loudness.loudnessTDep(startSkipIdx:end, :));
                            case '5PcEx'
                                N = prctile(loudness.loudnessTDep(startSkipIdx:end, :), 95);
                        end  % end of switch for aggregate
                    end  % end of if branch for aggregate
                end  % end of if branch for binaural
            case 'ISO 532-1'
                loudness = acousticLoudness(signal, sampleRateIn, 1,...
                                            'Method', 'ISO 532-1',...
                                            'TimeVarying', true,...
                                            'SoundField', soundFieldKey);
                switch aggregate
                    case 'PowAvg'
                        N = mean(loudness(startSkipIdx:end, :).^(1/log10(2)), 1).^log10(2);
                    case 'Mean'
                        N = mean(loudness(startSkipIdx:end, :), 1);
                    case 'Max'
                        N = max(loudness(startSkipIdx:end, :), [], 1);
                    case 'Min'
                        N = min(loudness(startSkipIdx:end, :), [], 1);
                    case '5PcEx'
                        N = prctile(loudness(startSkipIdx:end, :), 95, 1);
                end  % end of switch for aggregate
        end  % end of switch for method
    end


end

