function detectability = auralDetectability(signalTarget, sampleRateTarget, signalMasker, sampleRateMasker, axisTarget, axisMasker, timeStep, timeSkip, nOct, freqRange, outPlot)
% detectability = auralDetectability(signalTarget, sampleRateTarget,
%                                        signalMasker, sampleRateMasker,
%                                        axisTarget, axisMasker,
%                                        timeStep, timeSkip, nOct,
%                                        freqRange, outPlot)
%
% Returns aural detectability and discounted sound levels from input target
% source and masker signals based on the detectability model originally
% developed by Bolt, Beranek and Newman consulting engineers, and developed
% further by NASA, with discounting of target source sound levels using the
% technique developed by NASA (see References).
%
% Inputs
% ------
% signalTarget : vector or 2D matrix
%   Input target signal as single mono or stereo audio (sound pressure) 
%   signals.
%
% sampleRateTarget : integer
%   Sample rate (frequency) of the input target signal(s).
%
% signalMasker : vector or 2D matrix
%   Input masker signal(s) as single mono or stereo audio (sound pressure)
%   signals.
%
% sampleRateMasker : integer
%   Sample rate (frequency) of the input masker signal(s).
%
% axisTarget : integer (1 or 2, default: 1)
%   Time axis for the target signal(s) along which to determine detection.
%
% axisMasker : integer (1 or 2, default: 1)
%   Time axis for the masker signal(s) along which to determine detection.
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
%   Note that fractional-octave bands other than 1/3-octave have not been
%   validated.
%
% freqRange : vector (default: [20, 20000])
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
% dBSpecTarget : matrix
%                sound pressure level spectrogram for the input target
%                signal, with dimensions [timeOut, freqBands, targetChans]
%
% dBSpecMasker : matrix
%                sound pressure level spectrogram for the input masker
%                signal, with dimensions [timeOut, freqBands, maskerChans] 
%
% dBSpecDiscTarget : matrix
%                    sound pressure level spectrogram for the input target
%                    signal detectability-discounted, with dimensions
%                    [timeOut, freqBands, targetChans]
%
% dBATDepTarget : matrix or vector
%                 time-dependent A-weighted sound pressure level for the
%                 input target signal, with dimensions [timeOut, targetChans]
%
% dBATDepTarget : matrix or vector
%                 time-dependent A-weighted sound pressure level for the
%                 input masker signal, with dimensions [timeOut, targetChans]
%
% dBATDepDiscTarget : matrix or vector
%                     time-dependent A-weighted sound pressure level for the
%                     input target signal detectability-discounted, with
%                     dimensions [timeOut, targetChans]
%
% dBATDepDiscount : vector
%                   time-dependent detectability discount values for the
%                   A-weighted sound pressure levels of target signal vs
%                   masker signal, with dimensions [timeOut, targetChans]
%
% LAETarget : vector
%             overall A-weighted sound exposure level for each input target
%             signal channel
%
% LAEMasker : vector
%             overall A-weighted sound exposure level for each input masker
%             signal channel
%
% LAEDiscTarget : vector
%                 overall detectability-discounted A-weighted sound
%                 exposure level for each input target signal channel
%
% LAeqTarget : vector
%              overall A-weighted time-averaged sound level for each input
%              target signal channel
%
% LAeqMasker : vector
%              overall A-weighted time-averaged sound level for each input
%              masker signal channel
%
% LAeqDiscTarget : vector
%                  overall detectability-discounted A-weighted
%                  time-averaged sound level for each input target signal
%                  channel
%
% dBADiscount : vector
%               overall detectability discount for A-weighted levels for
%               each input target signal channel
%
% detectabilitydB : matrix
%                   detectability spectrogram for the input target signal
%                   vs masker signal (i.e., 10log_10[d']), with dimensions
%                   [timeOut, freqBands, targetChans]
%
% detectTDepMaxdB : matrix or vector
%                   band-maximum time-dependent detectability dB, with
%                   dimensions [timeOut, targetChans]
%
% detectTDepIntdB : matrix or vector
%                   band-integrated time-dependent detectability dB, with
%                   dimensions [timeOut, targetChans]
%
% detectMaxdB : vector
%               overall maximum detectability, dB, for each input target
%               signal channel
%
% detectIntdB : vector
%               overall integrated detectability, dB, for each input target
%               signal channel
%
% detectMaxPcdB : structure
%                 contains percent-based aggregated metrics for maximum
%                 detectability, comprising:
% Ex50 : vector
%        50% exceeded (median) dB, for each input target signal channel
%
% detectIntPcdB : structure
%                 contains percent-based aggregated metrics for maximum
%                 detectability, comprising:
% Ex50 : vector
%        50% exceeded (median) dB, for each input target signal channel
%
% freqBands : vector
%             1/3-octave band centre-frequencies for input freqRange
%
% timeOut : vector
%           the window centre times for the spectrograms
%
% Assumptions
% -----------
% The input signals are calibrated to units of acoustic pressure in Pascals
% (Pa).
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
% Date last modified: 20/03/2026
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
        signalTarget (:, :) double {mustBeReal}
        sampleRateTarget (1, 1) double {mustBePositive, mustBeInteger}
        signalMasker (:, :) double {mustBeReal}
        sampleRateMasker (1, 1) double {mustBePositive, mustBeInteger}
        axisTarget (1, 1) {mustBeInteger, mustBeInRange(axisTarget, 1, 2)} = 1
        axisMasker (1, 1) {mustBeInteger, mustBeInRange(axisMasker, 1, 2)} = 1
        timeStep (1, 1) double {mustBePositive} = 0.5
        timeSkip (1, 2) double {mustBeReal, mustBeNonnegative} = [0, 0]
        nOct (1, 1) {mustBeInteger, mustBePositive} = 3
        freqRange (1, 2) double {mustBeInRange(freqRange, 19, 20000)} = [19, 20000]
        outPlot {mustBeNumericOrLogical} = false
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Input checks and resampling
% adjust 20 Hz fMin down to 19 Hz
if min(freqRange) == 20
    freqRange = [19, max(freqRange)];
end

% orient axes
if axisTarget == 2
    signalTarget = signalTarget.';
end

if axisMasker == 2
    signalMasker = signalMasker.';
end

% check input signal sampling frequencies are equal, otherwise resample to
% match higher rate
if sampleRateMasker > sampleRateTarget
    sampleRate = sampleRateMasker;
    up = sampleRate/gcd(sampleRate, sampleRateTarget);
    down = sampleRateTarget/gcd(sampleRate, sampleRateTarget);
    signalTarget = resample(signalTarget, up, down);
elseif sampleRateTarget > sampleRateMasker
    sampleRate = sampleRateTarget;
    up = sampleRate/gcd(sampleRate, sampleRateMasker);
    down = sampleRateMasker/gcd(sampleRate, sampleRateMasker);
    signalMasker = resample(signalMasker, up, down);
else
    sampleRate = sampleRateTarget;
end

% Check time skip
% Overall time (s)
T = size(signalTarget, 1)/sampleRate;

if sum(timeSkip) > T
    error("Error: total 'timeSkip' argument is larger than signal duration. Check inputs.")
end

%% Signal processing
% Calculate 1/3-octave Leq
% 6th order Butterworth filter bank
octFiltBank = octaveFilterBank('Bandwidth', '1/3 octave',...
                               'SampleRate', sampleRate,...
                               'FrequencyRange', freqRange,...
                               'FilterOrder', 6, ...
                               'OctaveRatioBase', 10);

signalFiltBankTarget = octFiltBank(signalTarget);
signalFiltBankMasker = octFiltBank(signalMasker);

[leqtSpecTarget, ~] = timeAveragedLevel(signalFiltBankTarget, sampleRate, timeStep);
[leqtSpecMasker, ~] = timeAveragedLevel(signalFiltBankMasker, sampleRate, timeStep);

detectability = auralDetectFromLeq(leqtSpecTarget, leqtSpecMasker, timeStep, timeSkip, nOct, freqRange, outPlot);

%% Functions

function [leq, t] = timeAveragedLevel(p, fs, blockTime)
    % leq = timeAveragedLevel(p, fs, blockTime)
    %
    % Return block-time-averaged sound level
    %
    signalLen = size(p, 1);
    blockLen = round(blockTime*fs);
    if blockLen > signalLen
        error("Error: the period of the averaging block is less than the signal duration.")
    end
    
    % Calculate number of time blocks and block indices
    nBlocks = floor(signalLen/blockLen);
    iBlockStart = 1:blockLen:nBlocks*blockLen;
    
    % Get block centre indices
    if mod(blockLen, 2)  % odd
        iBlockCentre = iBlockStart + floor(blockLen / 2);
    else  % even
        iBlockCentre = iBlockStart + blockLen / 2;
    end
    
    % Calculate rolling window mean of squared pressure
    p2_mean = movmean(p.^2, blockLen, 1, 'Endpoints', 'fill');
    
    % Extract block means and calculate Leq
    leq = 10*log10(p2_mean(iBlockCentre, :)./4e-10);

    % Output time, starting at 0
    t = ((0:blockLen:(nBlocks - 1)*blockLen)/fs).';
end  % end of timeAveragedLevel function

end % end of auralDetectability function