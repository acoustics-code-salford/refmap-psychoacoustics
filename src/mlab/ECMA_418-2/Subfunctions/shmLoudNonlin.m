function loudNonlin = shmLoudNonlin(pRMS)
% loudNonlin = shmLoudNonlin(pRMS)
%
% Returns the result of the compressive nonlinearity transformation from
% RMS sound pressure to (unthresholded) loudness according to
% ECMA-418-2:2025 (the Sottek Hearing Model), Section 5.1.8, Equation 23.
%
% Inputs
% ------
% pRMS : array
%   RMS sound pressure value(s) [Pa]
%
% Returns
% -------
% loudNonlin : array
%   the corresponding nonlinearly-transformed value(s), the same size as
%   pRMS
%
% Assumptions
% -----------
% pRMS is a scalar or row vector. Because the eight thresholds p_ti
% (Table 2) are broadcast against pRMS along a leading singleton
% dimension, a column-vector pRMS would be broadcast incorrectly; reshape
% to a row vector (or call element-by-element) if needed.
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
        pRMS double {mustBeReal, mustBeNonnegative}
    end

%% Define constants

% Section 5.1.8 Equation 23/24 ECMA-418-2:2025 [c_N], [alpha]
cal_N = 0.0211668;
cal_Nx = 1.00132;  % calibration adjustment factor (Footnote 8 ECMA-418-2:2025)
a = 1.5;

% Section 5.1.8 Table 2 ECMA-418-2:2025 [p_ti], [nu_i] (nu_0 = 1 prepended)
p_threshold = 2e-5*10.^((15:10:85)/20).';
v = [1, 0.6602, 0.0864, 0.6384, 0.0328, 0.4068, 0.2082, 0.3994, 0.6434];

%% Signal processing

% Section 5.1.8 Equation 23 ECMA-418-2:2025
loudNonlin = cal_N*cal_Nx*(pRMS/20e-6).*prod((1 + (pRMS./p_threshold).^a).^((diff(v)/a)'), 1);

% end of function