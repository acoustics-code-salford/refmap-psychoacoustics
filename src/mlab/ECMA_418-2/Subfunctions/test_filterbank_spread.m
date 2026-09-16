% test_filterbank_spread.m
%
% Standalone diagnostic using ONLY the shared, pre-existing repository
% functions (shmResample, shmPreProc, shmOutMidEarFilter,
% shmAuditoryFiltBank, shmSignalSegment, shmBasisLoudness) - none of the
% fluctuation-strength-specific code is touched here at all. Purpose: to
% measure how widely the auditory filterbank's output actually spreads
% across critical bands for the calibration tone, independent of
% anything in acousticSHMFluctuation.m. If this spread is already wide
% using only shared/already-validated code, that rules out the
% fluctuation-strength-specific candidate-selection/weighting logic as
% the source of the wide per-band spread seen in diagLog, and points
% instead at the auditory filter's own frequency response (or the
% interpretation of "60 dB SPL" used to generate the test signal) as the
% place to keep looking. If the spread is narrow here, that would argue
% the wide spread previously observed is being introduced somewhere in
% the fluctuation-strength-specific code after all, despite the many
% checks already done.
%
% Usage
% -----
%   test_filterbank_spread(sine_1kHz_4Hz_60dB.wav, sine_1kHz_4Hz_60dB.sr)
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
% Date created: 16/09/2026

function test_filterbank_spread(p, sampleRateIn)

arguments (Input)
    p (:, :) double {mustBeReal}
    sampleRateIn (1, 1) double {mustBePositive, mustBeInteger}
end

if size(p, 2) > 1
    p = p(:, 1);  % first channel only
end

sampleRate48k = 48e3;
dz = 0.5;
halfBark = 0.5:dz:26.5;
nBands = length(halfBark);
bandCentreFreqs = shmBark2Hz(halfBark);

% Use ROUGHNESS's block/hop parameters here (16384/4096) purely because
% this is a lighter-weight, already-validated segmentation the shared
% code is routinely exercised with; the auditory filterbank's own
% frequency response does not depend on segmentation block size, so this
% is a fair like-for-like test of the shared filter chain regardless.
blockSize = 16384;
overlap = 0.75;
hopSize = (1 - overlap)*blockSize;

if sampleRateIn ~= sampleRate48k
    [p_re, ~] = shmResample(p, sampleRateIn);
else
    p_re = p;
end

pn = shmPreProc(p_re, blockSize, hopSize, true, false);
pn_om = shmOutMidEarFilter(pn, 'freeFrontal');
pn_omz = shmAuditoryFiltBank(pn_om, false);

basisLoudnessAll = zeros(nBands, 1);
for zBand = nBands:-1:1
    [signalSegmented, ~] = shmSignalSegment(pn_omz(:, zBand), 1, blockSize, overlap, 1, true);
    [~, bandBasisLoudness, ~] = shmBasisLoudness(signalSegmented, bandCentreFreqs(zBand));
    % use a block comfortably past the transient (skip the first few,
    % take one near the end but not the very last, matching the same
    % "past-transient, not the final edge block" logic used earlier)
    nBlocks = size(bandBasisLoudness, 2);
    useBlock = max(1, nBlocks - 3);
    basisLoudnessAll(zBand) = bandBasisLoudness(useBlock);
end

T = table((1:nBands).', bandCentreFreqs(:), basisLoudnessAll, ...
          'VariableNames', {'zBand', 'freqHz', 'basisLoudness'});
disp(T)

fig = figure;
semilogx(bandCentreFreqs, basisLoudnessAll, 'o-');
xlabel('Band centre frequency, Hz');
ylabel('Basis loudness, sone_{HMS}/Bark_{HMS}');
title('Auditory filterbank output spread (shared, untouched code)');
grid on;

end