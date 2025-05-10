function [soundNLevel, fm, fn] = loudNweighting(spectrL, octN, fLim, axisF)
% soundNLevel = loudNweighting(spectr, octN, fLim, axisF, tStep)
%
% Returns quasi-loudness N-weighted sound levels based on ISO 226:2023 and
% ISO 389-7:2019 equal loudness contours. Accepts varying fractional octave
% band resolution and interpolates level-dependent weighting functions
% accordingly over the 20 Hz-16 kHz range.
%
% Inputs
% ------
% spectrL : vector or 2D matrix
%           the continguous input sound level spectrum or spectra for
%           processing.
%
% octN : integer (default: 3)
%        the octave fraction denominator (e.g., 3 indicates 1/3-octave
%        spectra)
%
% fLim : vector (default: [20, 16e3])
%        the frequency limits for the input spectra (these are
%        automatically matched to nearest exact band centre-frequencies for
%        the selected octN resolution
%
% axisF : integer (default: 2)
%         the frequency band axis for series of spectra
% 
% Returns
% -------
% soundNLevel : vector or 2D matrix
%               the N-weighted sound level spectrum or spectra, oriented as
%               per the input.
%
% fm : vector
%      the fractional octave band exact mid-frequencies
%
% fn : vector
%      the fractional octave band nominal mid-frequencies
%
% Assumptions
% -----------
% The input is a fractional octave band unweighted sound level spectrum or
% series of sound level spectra in dB re 2e-5 Pa.
%
% Requirements
% ------------
% Curve Fitting Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 23/04/2025
% Date last modified: 23/04/2025
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
        spectrL (:, :) double {mustBeReal}
        octN (1, 1) double {mustBePositive, mustBeInteger} = 3
        fLim (1, 2) double {mustBeInRange(fLim, 19, 16000)} = [20, 16000]
        axisF (1, 1) {mustBeInteger, mustBeInRange(axisF, 1, 2)} = 2
    end

%% Load path
addpath(genpath(fullfile("refmap-psychoacoustics", "src", "mlab")))

%% Input checks
if axisF == 1
    spectrL = spectrL.';
    axisTpose = true;
else
    axisTpose = false;
end


%% Define constants
% ISO 226:2023 alpha values
alpha = [0.635, 0.602, 0.569, 0.537, 0.509, 0.482, 0.456, 0.433, 0.412,...
         0.391, 0.373, 0.357, 0.343, 0.330, 0.320, 0.311, 0.303, 0.300,...
         0.295, 0.292, 0.290, 0.290, 0.289, 0.289, 0.289, 0.293, 0.303,...
         0.323, 0.354];

% estimated extension of alpha for 16 kHz band
alphaExt = [alpha, 0.387];

% ISO 226:2023 Lu values
Lu = [-31.5, -27.2, -23.1, -19.3, -16.1, -13.1, -10.4, -8.20, -6.30,...
      -4.60, -3.20, -2.10, -1.20, -0.50, 0.00, 0.40, 0.50, 0.00, -2.70,...
      -4.20, -1.20, 1.40, 2.30, 1.00, -2.30, -7.20, -11.2, -10.9, -3.50];

% estimated extension of Lu for 16 kHz band
LuExt = [Lu, -2];

% diffuse field narrowband noise hearing thresholds for 18-25 year-olds with
% normal hearing from ISO 389-7:2019 (20 Hz - 16 kHz)
% NOTE: watch out for non-standard frequencies in the standard Table 1!
hearThresholdsDF3897 = [78.1, 68.7, 59.5, 51.1, 44.0, 37.5, 31.5, 26.5, 22.1,...
                        17.9, 14.4, 11.4, 8.4, 5.8, 3.8, 2.1, 1.0, 0.8, 1.9,...
                        0.5, -1.5, -3.1, -4.0, -3.8, -1.8, 2.5, 6.8, 8.4, 14.4,...
                        43.7];

%% Weighting inputs
% fractional octave frequencies
[fm, ~, ~, fn] = noctf(fLim, octN);

% interpolate functions over all frequencies
[fm3, ~, ~, ~] = noctf([20 16e3], 3);

fm3I1 = find(min(fm) == fm3);
fm3I2 = find(max(fm) == fm3);

if octN == 3
    alphaUse = alphaExt(fm3I1:fm3I2);
    LuUse = LuExt(fm3I1:fm3I2);
    hTUse = hearThresholdsDF3897(fm3I1:fm3I2);
else
    alphaFit = fit(fm3.', alphaExt.', 'cubicspline');
    LuFit = fit(fm3.', LuExt.', 'cubicspline');
    hTFit = fit(fm3.', hearThresholdsDF3897, 'cubicspline');
    alphaUse = alphaFit(fm);
    LuUse = LuFit(fm);
    hTUse = hTFit(fm);
end

alpha1k = repmat(alphaUse(fm==1e3), 1, size(alphaUse, 2));
Lu1k = repmat(LuUse(fm==1e3), 1, size(LuUse, 2));
hT1k = repmat(hTUse(fm==1e3), 1, size(hTUse, 2));

%% Signal processing
% calculate the weighting value for each spectral level
weighting = loudLfISO226_2023(alpha1k, Lu1k, hT1k, spectrL) - loudLfISO226_2023(alphaUse, LuUse, hTUse, spectrL);
% apply weighting to input spectra
soundNLevel = spectrL + weighting;

% if transposed input, transpose output
if axisTpose
   soundNLevel = soundNLevel.';
end

end  % end of function

