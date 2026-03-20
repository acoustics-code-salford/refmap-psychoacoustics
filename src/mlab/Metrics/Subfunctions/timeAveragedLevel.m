function [leq, timeOut] = timeAveragedLevel(p, sampleRate, blockTime, axisN)
% [leq, timeOut] = timeAveragedLevel(p, sampleRate, axisN, blockLen)
%   Return block-time-averaged (equivalent continuous) levels for input
%   time-series data.
%
%
%
%% Arguments validation
    arguments
        p (:, :, :) double {mustBeReal}
        sampleRate (1, 1) double {mustBePositive}
        blockTime (1, 1) double {mustBePositive}
        axisN (1, 1) double {mustBeInteger, mustBeInRange(axisN, 1, 2)} = 1
    end

%% Check inputs
if axisN == 2
    p = p.';
    % axisN = 1;
end

signalLen = size(p, 1);
blockLen = round(blockTime*sampleRate);
if blockLen > signalLen
    error("Error: the period of the averaging block is less than the signal duration.")
end

%% Processing

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
leq = 10*log10(p2_mean(iBlockCentre, :, :)./4e-10);

% Output time, starting at 0
timeOut = ((0:blockLen:(nBlocks - 1)*blockLen)/sampleRate).';
