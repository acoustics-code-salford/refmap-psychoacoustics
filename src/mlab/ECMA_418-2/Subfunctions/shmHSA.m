function [pHat, Elz, diagInfo] = shmHSA(fc, spectrumE, blockSize, sampleRate, nZerosStart, nZerosEnd, epsilon)
% [pHat, Elz, diagInfo] = shmHSA(fc, spectrumE, blockSize, sampleRate, nZerosStart, nZerosEnd, epsilon)
%
% Returns the High-resolution Spectral Analysis (HSA) estimate of the
% constant (DC) component and the complex spectral line amplitude(s) at
% the candidate modulation rate(s) fc, together with the corresponding
% HSA error function value, according to ECMA-418-2:2025 (the Sottek
% Hearing Model), Section 9.1.4.
%
% Inputs
% ------
% fc : row vector
%   candidate modulation rate(s) [Hz] of the Mc non-zero spectral line
%   pairs under consideration, fc = (fc_1, ..., fc_Mc). The zero
%   (constant) component is handled internally and must NOT be included
%   in fc. Negative values are accepted (the error function E_l,z is an
%   even function of each fc, since W+ is even and W- is odd in fc), so
%   that the Newton iteration of Section 9.1.7 may pass through zero and
%   the resulting f_c,1,opt < 0.125 Hz is then discarded by the caller as
%   the standard specifies, rather than raising an error.
%
% spectrumE : column vector
%   the s~b-point complex DFT spectrum P_E,l,z(k) of the windowed,
%   downsampled envelope (Section 9.1.4 Equation 121), full length
%   (k = 0, ..., s~b - 1), indexed with MATLAB index 1 corresponding to
%   k = 0
%
% blockSize : double
%   downsampled analysis block size s~b (Section 9.1.2)
%
% sampleRate : double
%   downsampled analysis sample rate r~s (Section 9.1.2)
%
% nZerosStart : double
%   number of zeros at the start of the envelope analysis window, n_zb,l,z
%
% nZerosEnd : double
%   number of zeros at the end of the envelope analysis window, n_ze,l,z
%
% epsilon : double (optional, default = 1e-12)
%   small constant substituted for the standard's epsilon_0 (see Note
%   below and in shmHSAWindowResponse.m)
%
% Returns
% -------
% pHat : column vector, length Mc + 1
%   the HSA-estimated complex spectral amplitudes: pHat(1) = phat_0,l,z
%   (real-valued constant part), pHat(2:end) = phat_fc,m,l,z (complex),
%   m = 1, ..., Mc, in the same order as the input fc. The spectral
%   line amplitudes are TWO-SIDED line amplitudes, i.e. an envelope
%   component a*cos(2*pi*fc*t + phi) is returned as
%   phat_fc = (a/2)*exp(1i*phi) (see the Note below)
%
% Elz : double
%   the HSA error function value E_l,z(fc) (Section 9.1.4 Equation 135)
%
% diagInfo : structure
%   diagnostic information about the linear solve, for use in isolating
%   numerical issues (this is not part of the standard's algorithm - it
%   is provided purely to support debugging/validation). Fields:
%     Mc       : number of non-zero candidate lines fitted (numel(fc))
%     KL       : the number of DFT bins used in the fit (Equation 125)
%     nUnknown : 2*Mc + 1, the number of unknowns in Formula (130)
%     rcondA   : rcond(A), the reciprocal condition number estimate of
%                the matrix in Formula (130). Values close to 0 indicate
%                a near-singular (poorly conditioned or rank-deficient)
%                system; this becomes structurally more likely as Mc
%                grows, since K_L is capped at 49 (Equation 125)
%                regardless of Mc, while the system has 2*Mc + 1
%                unknowns - so Mc approaching 24-25 can leave the system
%                exactly or nearly rank-deficient (49 equations against
%                up to 51 unknowns).
%     usedPinv : true if rcondA was low enough to trigger the pinv
%                fallback described in the Note below, in which case the
%                returned pHat is a minimum-norm least-squares solution
%                rather than a well-determined one
%
% Assumptions
% -----------
% fc contains only the non-zero candidate modulation rate(s); the
% constant (DC) part is always included as an additional unknown and
% must not be passed explicitly.
%
% Note
% ----
% This function implements the general Mc-line case (Section 9.1.4,
% Equations 121-135) directly, using MATLAB's backslash operator to
% solve the (2*Mc + 1)-by-(2*Mc + 1) linear system of Equation 130. This
% is also used for the single spectral line pair case (Mc = 1), for
% which Section 9.1.4.1 additionally provides a closed-form solution via
% Cramer's rule (Equations 136-142). The closed-form solution is
% mathematically identical to the general solution for Mc = 1 (it is
% simply an expanded, manual method of solving the same 3-by-3 instance
% of Formula (130)) so solving the general system with backslash is not
% an approximation or a shortcut: it produces the same result (to
% floating-point rounding) while keeping a single, simpler, and more
% easily verified code path for all values of Mc.
%
% Note also that Formula (126) defines the column of W corresponding to
% the "minus" component (index i = 2*m + 1) as the COMPLEX CONJUGATE of
% W-_E,l,z,fc,m(k) (Equation 129); this conjugation is essential to
% obtaining correct results and is easy to miss when transcribing the
% formulae, since Equation 129 itself does not show the conjugation
% (it is introduced only in Formula (126), and repeated in the a13/a23/b3
% terms of Formulae (138)-(139) for the Mc = 1 case).
%
% Note on the amplitude convention of phat_fc,m (Equation 123): the
% real-valued normal equations of Formulae (130)-(134) are exactly the
% least-squares fit of the model
%   Phat(k) = x_1*W_0(k) + sum_m [x_2m*W+_m(k) + j*x_2m+1*W-_m(k)]
% to P_E,l,z(k), which is the DFT of the windowed envelope model
%   x_1 + sum_m 2*Re((x_2m + j*x_2m+1)*exp(j*2*pi*fc,m*ntilde/rs~)).
% Hence (x_2m + j*x_2m+1) is the two-sided line amplitude (the line at
% +fc,m, with its conjugate at -fc,m), and 2*(x_2m + j*x_2m+1) would be
% the one-sided (cosine) amplitude. Equation 123 as printed
% (x_i = Re(phat)/2, Im(phat)/2) implies the one-sided convention, but
% Equations 159-160 (phat_0^2 + 2*sum(A_i) as the mean-square power of
% the harmonic complex, and sqrt(0.5*(...)) as the RMS sound pressure
% passed to the nonlinearity of Equation 23) and footnote 46 (the DFT
% line values of Equation 66 equal the HSA line values multiplied by the
% DFT length s~b) are only consistent with the two-sided convention.
% The two-sided convention is therefore used here, and has been
% verified against the reference implementation (HEAD acoustics
% ArtemiS v17) - see the note at the point of use below.
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
%
% Checked by:
% Date last checked:
%
%% Arguments validation
    arguments (Input)
        fc (1, :) double {mustBeReal, mustBeNonzero}
        spectrumE (:, 1) double
        blockSize (1, 1) double {mustBePositive}
        sampleRate (1, 1) double {mustBePositive}
        nZerosStart (1, 1) double {mustBeReal, mustBeNonnegative}
        nZerosEnd (1, 1) double {mustBeReal, mustBeNonnegative}
        epsilon (1, 1) double {mustBePositive} = 1e-12
    end

%% Define constants

Mc = numel(fc);  % number of non-zero candidate spectral line pairs
nCols = 2*Mc + 1;  % number of unknowns/columns, size of x [2Mc + 1]

% Section 9.1.4 Equation 125 [Delta f] and [K_L]
deltaF = sampleRate/blockSize;
KL = min(max(17, round(max(abs(fc))/deltaF) + 8), 49);

kIndices = (0:KL - 1).';  % zero-based DFT bin indices used in the fit

%% Signal processing

% Section 9.1.4 Equation 126 - build the matrix of window response
% column vectors W = (W_1, ..., W_2Mc+1)
W = zeros(KL, nCols);

% constant (DC) part, i = 1 [W_E,l,z,0(k)]
W(:, 1) = shmHSAWindowResponse(kIndices, 0, blockSize, sampleRate,...
                               nZerosStart, nZerosEnd, epsilon);

for mLine = 1:Mc
    windowPos = shmHSAWindowResponse(kIndices, fc(mLine), blockSize,...
                                     sampleRate, nZerosStart, nZerosEnd, epsilon);
    windowNeg = shmHSAWindowResponse(kIndices, -fc(mLine), blockSize,...
                                     sampleRate, nZerosStart, nZerosEnd, epsilon);

    % Equation 128 [W+_E,l,z,fc,m(k)], i = 2m
    W(:, 2*mLine) = windowPos + windowNeg;

    % Equation 129 [W-_E,l,z,fc,m(k)], conjugated per Equation 126, i = 2m+1
    W(:, 2*mLine + 1) = conj(windowPos - windowNeg);
end

% Section 9.1.4 Equations 132-133 - index sets I_R (and, by the same
% definition, J_R) identifying which columns are treated as the 'real'
% (cosine-type) part of the system
isRealCol = false(1, nCols);
isRealCol([1, 2:2:nCols]) = true;  % i = 1 or mod(i, 2) = 0
sameBucket = isRealCol.' == isRealCol;  % (i,j) both in, or both out of, I_R

% Section 9.1.4 Equation 131 - symmetric matrix A
% (the formula for a_ij is symmetric under exchange of i and j in both
% branches, so the full matrix can be built directly without separately
% enforcing symmetry)
Wr = real(W);
Wi = imag(W);
A = (Wr.'*Wr + Wi.'*Wi).*sameBucket + (Wi.'*Wr + Wr.'*Wi).*(~sameBucket);

% Section 9.1.4 Equation 134 - vector b
spectrumEk = spectrumE(kIndices + 1);  % PE,l,z(k) at the KL bins used
PEr = real(spectrumEk);
PEi = imag(spectrumEk);

b = zeros(nCols, 1);
b(isRealCol) = Wr(:, isRealCol).'*PEr + Wi(:, isRealCol).'*PEi;
b(~isRealCol) = Wr(:, ~isRealCol).'*PEi + Wi(:, ~isRealCol).'*PEr;

% Section 9.1.4 Equation 130 - solve A*x = b
% (a numerical robustness safeguard, not specified by the standard: fall
% back to the minimum-norm least-squares solution if A is close to
% singular, which can occur for pathological candidate frequency sets)
rcondA = rcond(A);
usedPinv = rcondA < 1e3*eps(class(A));
if usedPinv
    x = pinv(A)*b;
else
    x = A\b;
end

% Section 9.1.4 Equation 135 - error function
% (a_ii*x_i^2 summed with 2*sum_{i<j}(a_ij*xi*xj) is exactly x'*A*x for
% symmetric A, giving a compact, fully vectorised equivalent of Formula
% (135))
Elz = sum(abs(spectrumEk).^2) + x.'*A*x - 2*b.'*x;

% Section 9.1.4 Equation 123 (inverted) - recover the complex spectral
% amplitudes from the solution vector x.
%
% NOTE on the amplitude convention (see also the Note in the function
% help): the spectral line amplitude is returned as the TWO-SIDED line
% amplitude phat_fc,m = x_2m + j*x_2m+1 (i.e. the amplitude of the
% positive-modulation-rate line, whose conjugate sits at -fc,m), and NOT
% as 2*(x_2m + j*x_2m+1) as a literal inversion of Equation 123 would
% give. The two-sided convention is the one required by Equations 159
% and 160 (where phat_0^2 + 2*sum(A_i) is the mean-square of the
% envelope, so that sqrt(0.5*(...)) is the RMS of the band-pass signal
% fed to the nonlinearity of Equation 23) and by footnote 46 (the DFT
% results of Equation 66 equal the HSA results multiplied by the DFT
% length). Using the literal Equation 123 factor of 2 overestimates
% |phat|^2 by a factor of 4 and the harmonic-complex power by up to a
% factor of 2, and was found to reproduce a modulation-depth-dependent
% overestimation of fluctuation strength (approx. 1.8x at m = 1 rising
% to 3x at m = 0.25) relative to the reference (HEAD acoustics ArtemiS)
% results; with the two-sided convention the reference calibration
% signal yields F = 0.994 vacil_HMS (ArtemiS: 1.003) and per-band
% specific fluctuation strength agrees with ArtemiS to within a few %
% across 50-70 dB and modulation depths 25-100 %.
pHat = zeros(Mc + 1, 1);
pHat(1) = x(1);  % [phat_0,l,z], real-valued
for mLine = 1:Mc
    pHat(mLine + 1) = x(2*mLine) + 1i*x(2*mLine + 1);  % [phat_fc,m,l,z], two-sided line amplitude
end

% diagnostic information (not part of the standard - see Returns above)
diagInfo.Mc = Mc;
diagInfo.KL = KL;
diagInfo.nUnknown = nCols;
diagInfo.rcondA = rcondA;
diagInfo.usedPinv = usedPinv;

% end of function