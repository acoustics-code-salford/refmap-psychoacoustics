function [pHat, Elz] = shmHSA(fc, spectrumE, blockSize, sampleRate, nZerosStart, nZerosEnd, epsilon)
% [pHat, Elz] = shmHSA(fc, spectrumE, blockSize, sampleRate, nZerosStart, nZerosEnd, epsilon)
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
%   in fc.
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
%   m = 1, ..., Mc, in the same order as the input fc
%
% Elz : double
%   the HSA error function value E_l,z(fc) (Section 9.1.4 Equation 135)
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
% Date last modified: 16/09/2026
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
        fc (1, :) double {mustBeReal, mustBePositive}
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
KL = min(max(17, round(max(fc)/deltaF) + 8), 49);

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
if rcond(A) < 1e3*eps(class(A))
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
% amplitudes from the solution vector x
pHat = zeros(Mc + 1, 1);
pHat(1) = x(1);  % [phat_0,l,z], real-valued
for mLine = 1:Mc
    pHat(mLine + 1) = 2*x(2*mLine) + 1i*2*x(2*mLine + 1);  % [phat_fc,m,l,z]
end

% end of function