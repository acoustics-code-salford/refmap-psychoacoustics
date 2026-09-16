function windowResponse = shmHSAWindowResponse(kIndices, modRate, blockSize, sampleRate, nZerosStart, nZerosEnd, epsilon)
% windowResponse = shmHSAWindowResponse(kIndices, modRate, blockSize, sampleRate, nZerosStart, nZerosEnd, epsilon)
%
% Returns the analysis-window frequency response used by the High-
% resolution Spectral Analysis (HSA), according to ECMA-418-2:2025 (the
% Sottek Hearing Model), Section 9.1.4, Equation 127.
%
% Inputs
% ------
% kIndices : column vector
%   the (zero-based) DFT bin indices k at which to evaluate the window
%   response, k = 0, ..., K_L - 1
%
% modRate : double
%   the candidate modulation rate f_c,m [Hz] at which the window is
%   centred (may be zero, positive or negative)
%
% blockSize : double
%   the downsampled analysis block size s~b (Section 9.1.2)
%
% sampleRate : double
%   the downsampled analysis sample rate r~s (Section 9.1.2)
%
% nZerosStart : double
%   number of zeros at the start of the envelope analysis window, n_zb,l,z
%
% nZerosEnd : double
%   number of zeros at the end of the envelope analysis window, n_ze,l,z
%
% epsilon : double
%   small constant substituted for the standard's epsilon_0 (defined as
%   the smallest positive double such that 1 + epsilon_0 > 1), added to
%   avoid division by zero. See the Note in shmHSA.m for why the exact
%   value used here has no measurable effect on the result, and why the
%   repository-standard value of 1e-12 is used for consistency with the
%   other functions in this repository, in place of the true machine
%   epsilon.
%
% Returns
% -------
% windowResponse : column vector
%   the complex window response W_E,l,z,fc,m(k), same length as kIndices
%
% Assumptions
% -----------
% kIndices is a column vector of zero-based DFT bin indices, consistent
% with the zero-based indexing convention used throughout ECMA-418-2:2025.
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
        kIndices (:, 1) double {mustBeReal, mustBeNonnegative}
        modRate (1, 1) double {mustBeReal}
        blockSize (1, 1) double {mustBePositive}
        sampleRate (1, 1) double {mustBePositive}
        nZerosStart (1, 1) double {mustBeReal, mustBeNonnegative}
        nZerosEnd (1, 1) double {mustBeReal, mustBeNonnegative}
        epsilon (1, 1) double {mustBePositive} = 1e-12
    end

%% Signal processing

% number of samples with unity weight in the analysis window [n_active]
nActive = blockSize - nZerosEnd - nZerosStart;

% Section 9.1.4 Equation 127 [f_n(k)]
% (epsilon substitutes for the standard's epsilon_0, avoiding 0/0 when a
% candidate frequency falls exactly on a DFT bin centre - see Note above)
freqNorm = kIndices./blockSize - modRate./sampleRate + epsilon;

% NOTE: Equation 127 shows exp(-1i*2*pi*freqNorm...which seems to be
% incorrect
% As-written: windowResponse = exp(-1i*2*pi*freqNorm.*(blockSize - nZerosEnd + nZerosStart - 1))...
                 % .*sin(pi*freqNorm.*nActive)./sin(pi*freqNorm);
windowResponse = exp(-1i*pi*freqNorm.*(blockSize - nZerosEnd + nZerosStart - 1))...
                 .*sin(pi*freqNorm.*nActive)./sin(pi*freqNorm);


% end of function