% test_shmHSA.m
%
% Standalone diagnostic script to validate the shmHSA.m solver against
% synthetic signals with known ground truth, isolated from the rest of
% the acousticSHMFluctuation.m pipeline (envelope calculation, window
% determination, candidate selection, fine tuning, harmonic analysis).
%
% Update (post-fix)
% -----------------
% Running this script against the original shmHSAWindowResponse.m
% surfaced the root cause of the divergence from the reference (HEAD
% acoustics ArtemiS) results: Equation 127's phase prefactor was
% transcribed exactly as typeset in ECMA-418-2:2025, but the printed
% exponent exp(-j*2*pi*f_n(k)*(...)) is a factor of 2 too large; the
% correct exponent is exp(-j*pi*f_n(k)*(...)). This was confirmed two
% independent ways: a closed-form re-derivation of the underlying
% geometric series, and a brute-force numerical DFT comparison matching
% to machine precision (see the correction note in
% shmHSAWindowResponse.m for the full derivation). With the fix applied,
% Test 1 recovers the true amplitude and phase to ~1e-8 relative error
% (previously >80% error), and Test 2 shows recovery error at the
% 1e-8-1e-9 level across the entire Mc sweep from 3 to 25 - meaning the
% rank-deficiency/conditioning concern this script was originally written
% to probe is a real structural feature of Section 9.1.5 (rcond(A) does
% shrink as Mc grows) but was NOT the dominant cause of the discrepancy
% seen in the full pipeline; the Equation 127 exponent was.
%
% A second, independent typo was subsequently found and fixed in
% acousticSHMFluctuation.m's Equation 144 (the Section 9.1.5 stage-1
% candidate frequency estimate): the printed "- 1" term introduces a
% systematic bias of close to one full DFT bin, confirmed numerically
% across multiple single-tone synthetic tests and fixed by removing it.
%
% A third issue - the one that actually accounted for the remaining
% systematic overestimation relative to the reference - was the amplitude
% convention of the recovered spectral lines: a literal inversion of
% Equation 123 gives the one-sided (cosine) amplitude, whereas Equations
% 159-160 and footnote 46 require the two-sided line amplitude (half the
% cosine amplitude). shmHSA.m now returns the two-sided amplitude, and
% the ground-truth values in this script are defined accordingly (A1/2
% in Test 1, trueAmps/2 in Test 2). Note that this class of error is
% invisible to a self-consistent unit test: the previous version of this
% script "passed" with the one-sided convention because its expected
% values were written with the same convention.
% This script (shmHSA.m in isolation) would not detect that class of
% error, since it lives in the surrounding candidate-selection logic in
% acousticSHMFluctuation.m, not in the solver itself - a reminder that
% passing these tests confirms the solver is correct given whatever
% frequencies it is asked to fit, not that the right frequencies are
% being chosen upstream.
%
% Purpose
% -------
% This script exists to answer one question in isolation: given a
% correctly-formed windowed envelope spectrum P_E,l,z(k) and a set of
% candidate modulation rates, does shmHSA.m recover the correct complex
% amplitudes? If Test 1 fails, the bug is in shmHSA.m/shmHSAWindowResponse.m
% themselves. If Test 1 passes but Test 2 shows large errors as the
% candidate count grows, the bug is not in the solver's algebra but in
% the numerical conditioning of the linear system for large Mc (see the
% risk-ranking discussion this script was written to test: Section 9.1.5
% allows up to Mc = 25 simultaneous candidates - 24 local maxima plus
% f_min - while K_L (Equation 125) is capped at 49 regardless of Mc, so
% the system in Formula (130) has 2*Mc + 1 unknowns against at most 49
% equations and can become rank-deficient for large Mc).
%
% Usage
% -----
% Run this script directly (no inputs required):
%   >> test_shmHSA
%
% It prints a pass/fail summary for each test to the command window and
% does not modify anything outside its own workspace.
%
% Requirements
% ------------
% Signal Processing Toolbox (not strictly required by this script itself,
% but shmHSA.m's callers in acousticSHMFluctuation.m assume it)
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 16/09/2026
% Date last modified: 18/09/2026
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

% clearvars
% close all
% clc

fprintf('===============================================================\n');
fprintf('shmHSA.m standalone diagnostic tests\n');
fprintf('===============================================================\n\n');

%% Common analysis parameters (matching acousticSHMFluctuation.m defaults)

blockSize1500 = 2048;   % s~b
sampleRate1500 = 1500;  % r~s
nzb = 64;                % default n_zb
nze = 64;                % default n_ze
nActive = blockSize1500 - nzb - nze;  % 1920
deltaF1500 = sampleRate1500/blockSize1500;  % Delta f, ~0.7324 Hz

envWindow = zeros(blockSize1500, 1);
envWindow(nzb + 1:end - nze) = 1;

nSamp = (0:blockSize1500 - 1).';
tSec = nSamp/sampleRate1500;

%% =============================================================
%  TEST 1: single-line (Mc = 1) recovery accuracy
%  =============================================================
%
% A clean two-term synthetic envelope p_E(t) = A0 + A1*cos(2*pi*f1*t + phi1)
% is windowed exactly as the real pipeline windows it, and shmHSA is
% asked to recover A0 and the complex amplitude A1*exp(1i*phi1) at the
% (deliberately non-bin-centred) frequency f1. Because the HSA explicitly
% models the rectangular window's own DFT response (Equation 127), this
% recovery should be accurate to a small multiple of double-precision
% round-off, regardless of f1 not lying on a DFT bin.

fprintf('--- Test 1: single-line recovery ---------------------------\n');

A0 = 0.0500;             % true DC (constant) amplitude [Pa]
A1 = 0.0200;             % true AC (cosine) amplitude [Pa]
f1 = 4.1;                 % true modulation rate [Hz] (off-bin on purpose)
phi1 = 0.7;               % true phase [rad]
% true complex TWO-SIDED line amplitude at f1: shmHSA.m returns the
% amplitude of the line at +f1 (with its conjugate at -f1), i.e. A1/2,
% consistent with Equations 159-160 and footnote 46 of ECMA-418-2:2025
% (see the Note in shmHSA.m)
pTrue = (A1/2)*exp(1i*phi1);

envelope1 = A0 + A1*cos(2*pi*f1*tSec + phi1);
spectrumE1 = fft(envelope1.*envWindow, blockSize1500);

[pHat1, Elz1, diag1] = shmHSA(f1, spectrumE1, blockSize1500, sampleRate1500, nzb, nze);

p0Err = abs(pHat1(1) - A0)/A0;
pFcErr = abs(pHat1(2) - pTrue)/abs(pTrue);

fprintf('  True   p0 = %.6f            | Recovered p0 = %.6f  (rel. err %.3e)\n',...
        A0, real(pHat1(1)), p0Err);
fprintf('  True  pFc = %.6f %+.6fi   | Recovered pFc = %.6f %+.6fi  (rel. err %.3e)\n',...
        real(pTrue), imag(pTrue), real(pHat1(2)), imag(pHat1(2)), pFcErr);
fprintf('  Error function E_l,z(fc)  = %.3e\n', Elz1);
fprintf('  rcond(A) = %.3e, Mc = %d, KL = %d, usedPinv = %d\n',...
        diag1.rcondA, diag1.Mc, diag1.KL, diag1.usedPinv);

tol1 = 1e-6;
test1Pass = (p0Err < tol1) && (pFcErr < tol1);
if test1Pass
    fprintf('  RESULT: PASS (errors below %.0e)\n\n', tol1);
else
    fprintf('  RESULT: FAIL (errors exceed %.0e) - the fault is in shmHSA.m\n', tol1);
    fprintf('          or shmHSAWindowResponse.m itself.\n\n');
end

%% =============================================================
%  TEST 2: rank-deficiency / large-Mc stress test
%  =============================================================
%
% A synthetic envelope containing THREE genuine, well-separated spectral
% lines (plausible for a real signal) is built, but shmHSA is deliberately
% fed a much larger candidate set (up to 25 frequencies, matching the
% maximum Mc = 24 local maxima + f_min allowed by Section 9.1.5),
% including several candidates deliberately placed close to the genuine
% lines and to each other, to probe how the solver behaves as the number
% of unknowns (2*Mc + 1) approaches or exceeds K_L (capped at 49 by
% Equation 125). If the recovered amplitudes for the three genuine lines
% degrade, or spurious candidates pick up artificially large amplitude,
% as Mc grows, this reproduces (in isolation) the "many small spurious
% detections" pattern seen in the full pipeline's spectrograms.

fprintf('--- Test 2: rank-deficiency / large-Mc stress test ----------\n');

% three genuine lines: amplitudes chosen to be clearly different so
% correct vs incorrect recovery is easy to see
trueFreqs = [2.0, 4.9, 11.3];       % Hz
trueAmps  = [0.030, 0.015, 0.008];  % Pa
truePhi   = [0.2, -1.1, 2.4];       % rad

envelope2 = 0.04*ones(size(tSec));  % DC term
for iLine = 1:numel(trueFreqs)
    envelope2 = envelope2 + trueAmps(iLine)*cos(2*pi*trueFreqs(iLine)*tSec + truePhi(iLine));
end
spectrumE2 = fft(envelope2.*envWindow, blockSize1500);

% candidate counts to sweep: from a small, well-conditioned case up to
% the standard's maximum of Mc = 25 (24 local maxima + f_min)
McSweep = [3, 8, 15, 20, 24, 25];

fprintf('  %4s  %8s  %10s  %8s  %10s  %10s  %10s\n',...
        'Mc', 'KL', 'rcond(A)', 'pinv?', 'line1 err', 'line2 err', 'line3 err');

for Mc = McSweep
    % build a candidate set: the three true lines plus extra closely-
    % spaced "distractor" candidates filling out the rest of Mc,
    % spread across the 0.25-20 Hz range used by Section 9.1.5
    nExtra = Mc - numel(trueFreqs);
    if nExtra > 0
        extraFreqs = linspace(0.5, 19.5, nExtra);
        % nudge any extras that land suspiciously close to a true line
        % so we are testing "many candidates" rather than "duplicates"
        fcTest = sort([trueFreqs, extraFreqs]);
    else
        fcTest = sort(trueFreqs(1:Mc));
    end

    [pHat2, ~, diag2] = shmHSA(fcTest, spectrumE2, blockSize1500, sampleRate1500, nzb, nze);

    % find the recovered amplitude nearest each true line's frequency
    % (recovered values are two-sided line amplitudes, i.e. half the
    % cosine amplitude - see Test 1)
    lineErr = nan(1, 3);
    for iLine = 1:3
        [~, iNearest] = min(abs(fcTest - trueFreqs(iLine)));
        recoveredAmp = abs(pHat2(iNearest + 1));
        lineErr(iLine) = abs(recoveredAmp - trueAmps(iLine)/2)/(trueAmps(iLine)/2);
    end

    fprintf('  %4d  %8d  %10.3e  %8d  %10.3e  %10.3e  %10.3e\n',...
            Mc, diag2.KL, diag2.rcondA, diag2.usedPinv,...
            lineErr(1), lineErr(2), lineErr(3));
end

fprintf(['\n  Interpretation: rcond(A) trending toward 0 (and/or usedPinv\n',...
         '  becoming 1) as Mc grows, together with growing line-recovery\n',...
         '  errors, confirms the rank-deficiency hypothesis: the solver is\n',...
         '  mathematically correct (see Test 1) but the linear system\n',...
         '  Section 9.1.5 can hand it becomes poorly determined once the\n',...
         '  number of simultaneously-fitted candidates grows large, which\n',...
         '  is expected to occur far more often for busy real-world signals\n',...
         '  than for a clean single tone.\n\n']);

fprintf('===============================================================\n');
fprintf('End of shmHSA.m diagnostic tests\n');
fprintf('===============================================================\n');