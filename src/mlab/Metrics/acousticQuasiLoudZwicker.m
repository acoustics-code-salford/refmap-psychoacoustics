function loudness = acousticQuasiLoudZwicker(spectrL, fLim, axisF, soundField, adjustLoud)
% loudness = acousticQuasiLoudZwicker(spectrL, fLim, axisF, soundField,
%                                      adjustLoud)
%
% Returns quasi-loudness using spectral elements of ISO 532-1 Zwicker
% loudness model for arbitrary spectra. The temporal weighting of the
% model is disregarded. The input spectra must be 1/3 octave band
% resolution Leq, limited to the 25 Hz-12.5 kHz range.
%
% Optional modifications available comprise adjustments to the spectral
% levels to mimic the application of the ISO 532-3:2023 or
% ECMA-418-2:2025 outer-middle ear filter responses (in 1/3-octaves)
% instead of the ISO 532-1:2017 outer ear transmission, and
% approximated non-linear transformations more closely following
% ISO 532-3 or ECMA-418-2 instead of the ISO 532-1 original.
%
% Inputs
% ------
% spectrL : vector, 2D or 3D matrix
%   The contiguous input sound level 1/3-octave Leq spectrum or spectra for
%   processing. Must be orientated with time and frequency bands on
%   first two axes, and channels on third axis.
%
% fLim : vector (default: [25, 12.5e3])
%   The frequency limits for the input spectra (these are
%   automatically matched to nearest exact band centre-frequencies).
%
% axisF : integer optional (1 or 2, default: 2)
%   The frequency band axis for series of spectra.
%
% soundField : keyword string (default: 'freeFrontal')
%   Determines whether the 'freeFrontal' or 'diffuse' field stages
%   are applied in the outer-middle ear filter, or 'noOuter'
%   omits this filtering stage.
%
% adjustLoud : keyword string (default: 'iso5321')
%   Indicates whether to apply adjustments for the outer-middle ear filter
%   response and loudness transformations from ISO 532-3:2023 ('iso5323')
%   or ECMA-418-2:2025 ('ecma4182').
%   The default option ('iso5321') applies no adjustment, so follows ISO
%   532-1 more closely (for closer agreement with Zwicker's model).
%
% Returns
% -------
% loudness : structure
%   Contains the loudness output.
%
% loudness contains the following outputs:
%
% loudTDep :  vector or matrix
%   Time-dependent loudness arranged as [time(, channels)].
% 
% loudPowAvg : number or vector
%   Time-power-averaged loudness arranged as [loudness(, channels)].
%
% loud5pcEx : number or vector
%   95th percentile (5% exceeded) loudness
%   arranged as [loudness(, channels)].
%
% loudLevel :  vector or matrix
%   Time-dependent loudness level arranged as [time(, channels)].
%
% loudLvlPowAvg : number or vector
%   Time-power-averaged loudness level arranged as [loudness(, channels)].
%
% specLoudAvg : vector or matrix
%   Time-power-averaged specific loudness arranged as [bands(, channels)].
%
% specLoudAvg : matrix
%   Time-dependent specific loudness arranged as [time, bands(, channels)].
%
% barkAxis : vector
%   Critical band rates for specific loudness (0.1 dz intervals).
%
% freqInMid : vector
%   The 1/3-octave band exact mid-frequencies for the input spectra.
%
% freqInNom : vector
%   The 1/3-octave band nominal mid-frequencies for the input spectra.
%
%
% Assumptions
% -----------
% The input is a 1/3-octave band unweighted sound level spectrum or
% series of sound level spectra in dB re 2e-5 Pa.
%
% References
% ----------
%
% Zwicker loudness is defined in ISO 532-1:2017. Modifications are based on
% ISO 532-3:2023 and ECMA-418-2:2025.
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
% Date created: 23/04/2025
% Date last modified: 06/02/2026
% MATLAB version: 2023b
%
% Copyright statement: This file and code is part of work undertaken within
% the RefMap project (www.refmap.eu), and is subject to licence as detailed
% in the code repository
% (https://github.com/acoustics-code-salford/refmap-psychoacoustics)
%
% This file contains code adapted from the SQAT toolbox (itself being an
% adaptation of the AARAE code Loudness_ISO532_1.m by Ella Manor). SQAT is
% licensed as CC-NC-BY4 (https://creativecommons.org/licenses/by-nc/4.0/),
% while the AARAE code is BSD-3 (https://opensource.org/license/bsd-3-clause).
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
        spectrL (:, :, :) double {mustBeReal}
        fLim (1, 2) double {mustBeInRange(fLim, 25, 12600)} = [25, 12500]
        axisF (1, 1) {mustBeInteger, mustBeInRange(axisF, 1, 2)} = 2
        soundField (1, :) string {mustBeMember(soundField,...
                                               {'freeFrontal',...
                                                'diffuse',...
                                                'noOuter'})} = 'freeFrontal'
        adjustLoud (1, :) string {mustBeMember(adjustLoud,...
                                               {'iso5321',...
                                                'iso5323',...
                                                'ecma4182'})} = 'iso5321'
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Input checks
if ndims(spectrL) == 3
    inChans = size(spectrL, 3);
elseif ismatrix(spectrL)
    inChans = 1;
else
    error("Input spectrL must not have more than 3 dimensions.")
end

if axisF == 1
    if inChans == 1
        spectrL = spectrL.';
    else
        spectrL = permute(spectrL, [2, 1, 3]);
    end
end

% fractional octave frequencies
[freqInMid, ~, ~, freqInNom] = noctf(fLim, 3);
[fmAll, ~, ~, ~] = noctf([25, 12500], 3);

if length(freqInMid) ~= size(spectrL, 2)
    error("The frequency band range of the input spectra must correspond with the input band limits 'fLim'.")
end

%% Define constants
% ISO 532-1 Table A.3: low-frequency (25-250 Hz) weights
lowFWeights = [-32, -24, -16, -10, -5, 0, -7, -3, 0, -2, 0;
               -29, -22, -15, -10, -4, 0, -7, -2, 0, -2, 0;
               -27, -19, -14,  -9, -4, 0, -6, -2, 0, -2, 0;
               -25, -17, -12,  -9, -3, 0, -5, -2, 0, -2, 0;
               -23, -16, -11,  -7, -3, 0, -4, -1, 0, -1, 0;
               -20, -14, -10,  -6, -3, 0, -4, -1, 0, -1, 0;
               -18, -12,  -9,  -6, -2, 0, -3, -1, 0, -1, 0;
               -15, -10,  -8,  -4, -2, 0, -3, -1, 0, -1, 0];

% insert zeros to match all bands and tranpose to ease later calcs
lowFWeights = [lowFWeights, zeros([size(lowFWeights, 1), length(fmAll) - size(lowFWeights, 2)])].';

lowFLRange = [45; 55; 65; 71; 80; 90; 100; 120].';  % transposed for ease of calc

% ISO 532-1 Table A.4: (outer) ear transmission characteristic for free
% field incidence
earTmission  = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.5, -1.6,...
                -3.2, -5.4, -5.6, -4, -1.5, 2, 5, 12];

% ISO 532-1 Table A.5: diffuse field incidence adjustment
earDFAdjust = [0, 0, 0.5, 0.9, 1.2, 1.6, 2.3, 2.8, 3, 2,...
               0, -1.4, -2, -1.9, -1, 0.5, 3, 4, 4.3, 4];

% ISO 532-1 Table A.7: critical bandwidth vs. 1/3-octave bandwidth
% adjustment
critBWAdjust = [-0.25, -0.6, -0.8, -0.8, -0.5, 0, 0.5, 1.1, 1.5, 1.7,...
                1.8, 1.8, 1.7, 1.6, 1.4, 1.2, 0.8, 0.5, 0, -0.5];

switch adjustLoud
    case 'iso5321'
        % ISO 532-1 Table A.6: critical band level at threshold in quiet
        levelThresQ = [30, 18, 12, 8, 7, 6, 5, 4,...
                       3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3];

    case 'iso5323'
        % ISO 523-3 Table 2 (extended to 25 Hz)
        levelThresQTOB = [28.18, 28.18, 28.18, 28.18, 23.9, 19.2, 15.68,...
                          12.67, 10.09, 8.08, 6.3, 5.3, 4.5, 3.63, 3.63,...
                          3.63, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63,...
                          3.63, 3.63, 3.63, 3.63, 3.63];

        levelThresQCB = [10*log10(sum(10.^(levelThresQTOB(1:6)/10))),...
                         10*log10(sum(10.^(levelThresQTOB(7:9)/10))),...
                         10*log10(sum(10.^(levelThresQTOB(10:11)/10))),...
                         levelThresQTOB(12:end)];
    
        levelThresQ = levelThresQCB - critBWAdjust;

    case 'ecma4182'
        % ISO 523-3 Table 2 (extended to 25 Hz)    
        levelThresQTOB = [45, 39, 33, 28.18, 23.9, 19.2, 15.68,...
                          12.67, 10.09, 8.08, 6.3, 5.3, 4.5, 3.63, 3.63,...
                          3.63, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63,...
                          3.63, 3.63, 3.63, 3.63, 3.63];
    
        levelThresQCB = [10*log10(sum(10.^(levelThresQTOB(1:6)/10))),...
                         10*log10(sum(10.^(levelThresQTOB(7:9)/10))),...
                         10*log10(sum(10.^(levelThresQTOB(10:11)/10))),...
                         levelThresQTOB(12:end)];
    
        levelThresQ = levelThresQCB - critBWAdjust;

end

% ISO 532-1 Table A.8: critical band number (rate) - addition of 1e-4 is as
% as per Annex A.4
barkN = [0.9, 1.8, 2.8, 3.5, 4.4, 5.4, 6.6, 7.9, 9.2, 10.6, 12.3,...
         13.8, 15.2, 16.7, 18.1, 19.3, 20.6, 21.8, 22.7, 23.6, 24] + 1e-4;

% mapped to full frequency band range (omitting top band used only for
% slope)
barkNMap = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 1.8, 1.8, 1.8, 2.8, 2.8, 3.5,...
            4.4, 5.4, 6.6, 7.9, 9.2, 10.6, 12.3, 13.8, 15.2, 16.7, 18.1,...
            19.3, 20.6, 21.8, 22.7, 23.6] + 1e-4;

% ISO 532-1 Table A.9: specific loudness range bounds corresponding with
% loudness upper slope steepness
specNSlopeBound = [21.5; 18.0; 15.1; 11.5; 9.00; 6.10; 4.40; 3.10; 2.13;...
                   1.36; 0.82; 0.42; 0.00; 0.30; 0.22; 0.15; 0.00; 0.10;...
                   0.035; 0.00];

% ISO 532-1 Table A.9: specific loudness range bounds corresponding with
% loudness upper slope steepness
specNSlope = [13.0, 8.20, 6.30, 5.50, 5.50, 5.50, 5.50, 5.50;
              9.00, 7.50, 6.00, 5.10, 4.50, 4.50, 4.50, 4.50;
              7.80, 6.70, 5.60, 4.90, 4.40, 3.90, 3.90, 3.90;
              6.20, 5.40, 4.60, 4.00, 3.50, 3.20, 3.20, 3.20;
              4.50, 3.80, 3.60, 3.20, 2.90, 2.70, 2.70, 2.70;
              3.70, 3.00, 2.80, 2.35, 2.20, 2.20, 2.20, 2.20;
              2.90, 2.30, 2.10, 1.90, 1.80, 1.70, 1.70, 1.70;
              2.40, 1.70, 1.50, 1.35, 1.30, 1.30, 1.30, 1.30;
              1.95, 1.45, 1.30, 1.15, 1.10, 1.10, 1.10, 1.10;
              1.50, 1.20, 0.94, 0.86, 0.82, 0.82, 0.82, 0.82;
              0.72, 0.67, 0.64, 0.63, 0.62, 0.62, 0.62, 0.62;
              0.59, 0.53, 0.51, 0.50, 0.42, 0.42, 0.42, 0.42;
              0.40, 0.33, 0.26, 0.24, 0.24, 0.22, 0.22, 0.22;
              0.27, 0.21, 0.20, 0.18, 0.17, 0.17, 0.17, 0.17;
              0.16, 0.15, 0.14, 0.12, 0.11, 0.11, 0.11, 0.11;
              0.12, 0.11, 0.10, 0.08, 0.08, 0.08, 0.08, 0.08;
              0.09, 0.08, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05;
              0.06, 0.05, 0.03, 0.02, 0.02, 0.02, 0.02, 0.02];

specNSlope = [specNSlope, repelem(specNSlope(:, end), 1,...
                                  length(barkN) - size(specNSlope, 2))];

% % level difference between 2023 ISO 226 equal loudness contours (10-100
% % phon, 25-12500 Hz) and 1987 ISO 226 contours
% d2023v1987 = [ 8.4, 10.4, 11.0, 11.0, 10.6,  9.8,  8.5, 6.9, 4.8,  0.0;
%                  9.5, 12.0, 12.8, 12.8, 12.4, 11.5, 10.1, 8.4, 6.3,  3.6;
%                 10.2, 12.8, 13.7, 13.7, 13.1, 12.1, 10.6, 8.8, 6.6,  4.0;
%                 10.4, 13.1, 13.9, 13.8, 13.0, 11.9, 10.4, 8.5, 6.2,  3.5;
%                 10.3, 13.1, 13.8, 13.5, 12.7, 11.5,  9.9, 7.9, 5.6,  2.9;
%                 10.0, 12.9, 13.5, 13.3, 12.5, 11.3,  9.7, 7.7, 5.4,  2.7;
%                  9.5, 12.4, 13.2, 13.1, 12.3, 11.2,  9.7, 7.8, 5.5,  2.7;
%                  9.1, 12.0, 12.9, 12.9, 12.3, 11.2,  9.8, 7.9, 5.6,  2.9;
%                  8.2, 11.1, 12.2, 12.3, 12.0, 11.0,  9.7, 7.9, 5.7,  3.1;
%                  7.0,  9.9, 11.2, 11.6, 11.3, 10.5,  9.4, 7.8, 5.9,  3.3;
%                  5.9,  8.8, 10.0, 10.6, 10.5,  9.9,  9.0, 7.5, 5.7,  3.3;
%                  4.9,  7.5,  8.8,  9.4,  9.4,  9.1,  8.2, 6.9, 5.3,  3.2;
%                  3.5,  5.8,  7.1,  7.7,  7.8,  7.5,  6.8, 5.9, 4.5,  2.6;
%                  2.3,  4.2,  5.3,  5.8,  6.0,  5.9,  5.5, 4.7, 3.6,  2.2;
%                  1.2,  2.5,  3.3,  3.7,  3.9,  3.8,  3.4, 3.0, 2.2,  1.3;
%                  0.1,  0.7,  1.0,  1.3,  1.4,  1.4,  1.3, 1.1, 0.8,  0.4;
%                  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0;
%                  1.5,  1.5,  1.5,  1.6,  1.8,  2.0,  2.3, 2.7, 3.1,  3.7;
%                  1.8,  2.3,  2.6,  3.0,  3.3,  3.7,  4.2, 4.9, 5.7,  6.7;
%                  0.1,  0.5,  0.7,  1.1,  1.6,  2.3,  3.1, 4.0, 5.2,  6.5;
%                 -0.5,  0.1,  0.4,  0.9,  1.5,  2.3,  3.3, 4.4, 5.8,  7.3;
%                  0.4,  1.2,  1.7,  2.4,  3.0,  3.9,  5.0, 6.2, 7.5,  9.1;
%                  1.7,  2.8,  3.5,  4.2,  4.9,  5.9,  6.8, 8.0, 9.2, 10.7;
%                  2.8,  3.7,  4.5,  5.2,  5.8,  6.6,  7.4, 8.3, 9.2, 10.3;
%                  2.2,  3.0,  3.8,  4.3,  4.9,  5.5,  5.9, 6.5, 7.0,  7.4;
%                  0.5,  2.0,  3.1,  3.8,  4.4,  4.8,  5.0, 5.0, 4.9,  4.5;
%                  2.7,  5.5,  7.0,  7.9,  8.3,  8.2,  7.8, 7.0, 5.7,  0.0;
%                  6.9, 10.3, 12.0, 12.7, 12.5, 11.6,  9.7, 7.0, 3.2,  0.0];
% 
% phonRange = 10:10:100;

switch adjustLoud
    case 'iso5323'
        % Outer and middle ear filter functions corresponding with ISO 532-3,
        % in 1/3-octave bands (25 - 12500 Hz)
        % NOTE: these are determined from direct computation by the ISO
        % software, rather than using the tabulated values in Table 1
        outMidFree = [-31.23, -25.57, -21.36, -18.53, -16.02, -14.15,...
                      -12.55, -10.89, -9.355, -7.815, -6.481, -4.742,...
                      -3.240, -2.028, -0.855, -0.046, -0.157, -1.035,...
                       0.539,  3.593,  6.428,  8.058,  7.468,  3.887,...
                      -2.911, -9.555, -11.46, -11.18];
        
        outMidDiffuse = [-31.23, -25.57, -21.36, -18.53, -16.02, -14.14,...
                         -12.55, -10.89, -9.360, -7.941, -6.865, -5.089,...
                         -3.268, -2.070, -1.097,  0.268,  1.037,  0.965,...
                          1.170,  1.868,  4.764,  7.100,  6.142,  3.738,...
                         -0.882, -3.676, -5.503, -10.42];

        midOnly = [-31.22, -25.56, -21.36, -18.53, -16.02, -14.14, -12.56,...
                   -11.01, -9.652, -8.341, -7.402, -6.146, -4.829, -3.821,...
                   -3.314, -2.658, -2.777, -4.394, -5.968, -8.415, -9.763,...
                   -7.194, -6.663, -7.153, -9.931, -11.87, -10.84, -14.64];
        
        
    case 'ecma4182'
        % Outer and middle ear filter functions corresponding with ECMA-418-2:2024,
        % in 1/3-octave bands (25 - 12500 Hz)
        % NOTE: these are determined from direct computation using the
        % ECMA-defined filters
        outMidFree = [-22.29, -20.24, -18.31, -16.49, -14.68, -12.76,...
                      -11.14, -9.602, -8.242, -6.950, -5.702, -4.391,...
                      -2.898, -1.337,  0.180,  1.267,  0.033, -2.712,...
                      -1.747,  1.301,  3.552,  5.243,  6.166,  3.775,...
                      -2.485, -8.339, -10.08, -10.65];
        
        outMidDiffuse = [-22.29, -20.24, -18.31, -16.50, -14.68, -12.76,...
                         -11.15, -9.611, -8.257, -6.977, -5.751, -4.483,...
                         -3.090, -1.765, -0.768, -0.174,  0.232,  1.036,...
                          2.276,  3.569,  4.723,  5.871,  6.523,  3.994,...
                         -2.353, -8.266, -10.05, -10.63];
        
        midOnly = [-22.29, -20.24, -18.31, -16.50, -14.68, -12.76, -11.15,...
                   -9.611, -8.255, -6.971, -5.735, -4.445, -3.000, -1.572,...
                   -0.405,  0.369,  0.627,  0.421, -0.161, -1.024, -2.029,...
                   -3.051, -3.998, -4.768, -5.321, -5.664, -6.011, -7.418];
end

% Calibration constants (to ensure 1 sone corresponds with 1 kHz sinusoid at
% 40 dB in free-field)
switch adjustLoud
    case 'iso5321'
        calN = 1.192905697913430;
    case 'iso5323'
        calN = 1.263458171797368;
    case 'ecma4182'
        calN = 0.973779011817019;
end

% loudness transform threshold factor
threshFact = 0.25;

% loudness transform parameters (per critical band for ecma4182 option)
switch adjustLoud
    case 'ecma4182'
        tFormParams = [0.135021389	0.244761383, 0.20220873, 0.110702393,...
                       0.11654344, 0.112055262, 0.111001403, 0.122608561,...
                       0.109999749, 0.127299369, 0.116266959, 0.120431342,...
                       0.117903886, 0.10878395, 0.111512081, 0.114087954,...
                       0.116399176, 0.124184614, 0.131453002, 0.13591126;
                       8.923151571	0.29350401, 19.8856405, 9.761047979,...
                       6.104402159, 5.817516445, 23.65165966, 12.81430299,...
                       1.000021211, 24.26023815, 18.20581153, 15.14548704,...
                       9.634984752, 2.781608005, 8.063117329, 2.644819589,...
                       4.541834635, 10.46488604, 11.56702361, 11.10639455;
                       0.102031687, 0.10005271, 0.353280739, 0.819383619,...
                       0.846056284, 0.834337877, 0.640216654, 0.56709966,...
                       0.258048614, 0.66150985, 0.636190228, 0.717189801,...
                       0.834529929, 0.538910228, 0.890175983, 0.750378202,...
                       0.841244974, 0.812446455, 0.764877702, 0.790949093;
                       0.174802203, 0.143431283, 0.158014924, 0.192552183,...
                       0.191916773, 0.193541003, 0.195078403, 0.206931212,...
                       0.21013939, 0.200017825, 0.204760425, 0.233484617,...
                       0.19291123, 0.193407819, 0.187964204, 0.186507062,...
                       0.187365692, 0.187634862, 0.183894524, 0.172263166;
                      -0.065213305, -0.050093227, -0.033559524, -0.013707902,...
                      -0.014049748, -0.010297943, -0.010068442, -0.034871953,...
                      -0.022253104, -0.023742289, -0.022854561, -0.052904809,...
                      -0.009671036, -0.00533088, -0.002584325, -0.00367744,...
                      -0.007091664, -0.013850574, -0.016955836, -0.012684027;
                       82.27410377, 83.66590474, 84.15095885, 86.72797373,...
                       89.60190913, 88.67108652, 90.80936999, 125.2422638,...
                       119.4977957, 115.7697211, 116.2149252, 139.7546844,...
                       104.1361001, 95.73648912, 90.22242857, 93.35155915,...
                       99.94515072, 104.5910689, 103.9748482, 97.40956424;
                      -104.0115067, -47.08112111, -39.73853804, -24.79047747,...
                      -27.43525556, -24.27103556, -22.82474449, -60.99357938,...
                      -35.01262963, -39.10952994, -36.08846324, -35.76127856,...
                      -18.94203912, -12.42774782, -2.175130978, -2.438729452,...
                      -9.021924348, -15.31551128, -16.79954261, -12.36052383];


    case 'iso5323'
        tFormParams = [0.006, 0.034349337, 0.048173086, 0.046216707,...
                       0.056201012, 0.060524164, 0.064877972, 0.05129259,...
                       0.050999899, 0.056137975, 0.05912088, 0.059707345,...
                       0.061722822, 0.065657043, 0.091045287, 0.063191862,...
                       0.063480307, 0.073296674, 0.087537914, 0.095252625;
                       -0.25, -0.25, -0.25, -0.25, -0.3, -0.35, -0.25,...
                       -0.25, -0.25, -0.35, -0.304181474, -0.25, -0.3,...
                       -0.3, -0.3, -0.25, -0.35, -0.35, -0.3, -0.25;
                       0.5, 0.750000835, 0.748952577, 0.761347657,...
                       0.100035651, 0.100028378, 0.1, 0.1, 0.1, 0.100709475,...
                       0.096641999, 0.091373494, 0.100034744, 0.100030315,...
                       0.106857855, 0.09, 0.1, 0.08, 0.07, 0.1;
                       0.37, 0.302106926, 0.286776577, 0.288487659,...
                       0.279933692, 0.276340571, 0.273251194, 0.282216459,...
                       0.293978767, 0.275533377, 0.274355827, 0.272195311,...
                       0.269920312, 0.266103317, 0.241094842, 0.312964369,...
                       0.307714221, 0.25027212, 0.236054634, 0.275567116;
                       -0.012845595, -0.010511666, -0.00930198, -0.007484265,...
                       -0.011406429, -0.013873118, -0.0113972, -0.004562275,...
                       -0.021341126, -0.006241193, -0.00593841, -0.006148646,...
                       -0.003945805, -0.003384011, -0.021550635, -0.046389599,...
                       -4.62e-02, -1.04e-07, -1.05e-07, -0.054912583;
                       80.27522984, 97.40586961, 95.43366477, 94.77467445,...
                       95.64933715, 99.90313571, 97.44108654, 92.70782527,...
                       119.9999655, 94.06052415, 94.65498408, 91.95400275,...
                       88.5291194, 87.20458832, 65.4998204, 121.1877277,...
                       159.3429148, 165.7563846, 186.2492062, 167.8920953;
                       -41.3170001, -17.21325233, -15.17941803, -13.28052751,...
                       -18.4905788, -22.53931934, -20.0072265, -10.00473731,...
                       -43.0891087, -14.2439562, -14.98234127, -17.55498083,...
                       -12.5507699, -13.22727802, -58.98387763, -5.928168979,...
                       -6.138660497, -65.86515765, -68.87677816, -7.825806835];

end

%% Signal processing

% frequency band limiting indices
fmAllI1 = find(min(freqInMid) == fmAll);
fmAllI2 = find(max(freqInMid) == fmAll);
fmShift1 = fmAllI1 - 1;

% calculate number of critical bands for output
barkNMapOut = barkNMap(fmAllI1:fmAllI2);
barkNOut = unique(barkNMapOut);
barkCount = length(barkNOut);
barkI1 = find(min(barkNOut) == barkN);
barkI2 = find(max(barkNOut) == barkN);

[~, iBarkMin] = min(round(barkN, 1));
if iBarkMin == 1
    barkMin = 0.1;
else
    barkMin = barkN(iBarkMin - 1) + 0.1;
end
barkMax = max(round(barkN, 1));
barkAxis = barkMin:0.1:barkMax;
barkAxisN = length(barkAxis);

% number of time steps in spectral series
nTimeSteps = size(spectrL, 1);

% loop over channels
for chan = inChans:-1:1
    
    % low-frequency weighting and adjustment for equal loudness contour
    % differences
    levelWeighted = spectrL(:, :, chan);

    rangeInts = length(lowFLRange);
    % rangePhons = length(phonRange);
    
    % loop though 1/3-oct bands
    for kk = fmAllI1:fmAllI2
    
        % only proceed if any band weightings are non-zero
        if any(lowFWeights(kk, :)) && strcmp(adjustLoud, 'iso5321')
    
            % repeat band levels and apply all weights
            spectrLRepWt = repmat(spectrL(:, kk, chan), 1, rangeInts) + lowFWeights(kk, :);
    
            % check which weighted levels meet corresponding range criterion
            meetCriterion = spectrLRepWt <= lowFLRange;
    
            % get column indices for lowest criterion match
            [~, lowCriterionMatch] = max(meetCriterion == 1, [], 2);
    
            % repeat weights matrix for each time step
            lowFWeightsRep = repmat(lowFWeights(kk, :), nTimeSteps, 1);
    
            % weight levels according to level criteria
            levelWeighted(:, kk) = spectrL(:, kk, chan)...
                                   + lowFWeightsRep(sub2ind(size(lowFWeightsRep),...
                                   (1:nTimeSteps).', lowCriterionMatch));
        end  % end of if branch for non-zero weighting

        % % if adjusting for difference between 1987 vs 2023 equal loudness
        % % contours
        % if strcmp(adjustLoud, 'iso226')
        %     levelWeightedRep = repmat(levelWeighted(:, kk), 1, rangePhons);
        % 
        %     % check which levels meet corresponding criterion
        %     meetCriterion = levelWeightedRep <= phonRange;
        % 
        %     % get column indices for lowest criterion match
        %     [~, dCriterionMatch] = max(meetCriterion == 1, [], 2);
        % 
        %     % repeat differences matrix for each time step
        %     d2023v1987Rep = repmat(d2023v1987(kk, :), nTimeSteps, 1);
        % 
        %     % adjust levels according to matched criteria
        %     levelWeighted(:, kk) = levelWeighted(:, kk)...
        %                            - d2023v1987Rep(sub2ind(size(d2023v1987Rep),...
        %                            (1:nTimeSteps).', dCriterionMatch));
        % 
        % end  % end of if branch for 1987 vs 2023 equal loudness adjustment
    
    end  % end of for loop through bands
    
    % apply if using ECMA outer-middle ear filter function
    if strcmp(adjustLoud, 'ecma4182')
        switch soundField
            case "freeFrontal"
                levelWeighted = levelWeighted + outMidFree(fmAllI1:fmAllI2);
            case "diffuse"
                levelWeighted = levelWeighted + outMidDiffuse(fmAllI1:fmAllI2);
            case "noOuter"
                levelWeighted = levelWeighted + midOnly(fmAllI1:fmAllI2);
        end
    end

    % Calculate levels in critical bands
    upperBarkLen = length(levelWeighted(1, max(fmAllI1, 12):fmAllI2 - fmShift1));
    lvlWeightCB = [zeros([nTimeSteps, barkCount - upperBarkLen]), levelWeighted(:, max(fmAllI1, 12):fmAllI2 - fmShift1)];
    CB1 = false;
    CB2 = false;
    
    % if first critical band is included in frequency band range
    if any(barkNMapOut == barkN(1)) && fmAllI2 >= 6
    
        lvlWeightCB1 = 10*log10(sum(10.^(levelWeighted(:, max(fmAllI1,...
                                                              1):6)/10), 2));
        CB1 = true;
        lvlWeightCB(:, 1) = lvlWeightCB1;
    end
    
    % if second critical band is included in frequency band range
    if any(barkNMapOut==barkN(2)) && fmAllI2 >= 9
        
        lvlWeightCB2 = 10*log10(sum(10.^(levelWeighted(:, max(fmAllI1,...
                                                              7):9)/10), 2));
        CB2 = true;
    
        if CB1
            lvlWeightCB(:, 2) = lvlWeightCB2;
        else
            lvlWeightCB(:, 1) = lvlWeightCB2;
        end
    end
    
    % if third critical band is included in frequency band range
    if any(barkNMapOut==barkN(3)) && fmAllI2 >= 11
    
        lvlWeightCB3 = 10*log10(sum(10.^(levelWeighted(:, max(fmAllI1,10):11)/10), 2));
    
        if CB1
            lvlWeightCB(:, 3) = lvlWeightCB3;
        elseif CB2
            lvlWeightCB(:, 2) = lvlWeightCB3;
        else
            lvlWeightCB(:, 1) = lvlWeightCB3;
        end
    end
    
    % Loudness calculation
    % apply outer ear filter (if not using the ECMA outer-middle ear filter)
    if ~strcmp(adjustLoud, 'ecma4182')
        switch soundField
            case "freeFrontal"
                excitation = lvlWeightCB - earTmission(barkI1:barkI2);
            case "diffuse"
                excitation = lvlWeightCB - earTmission(barkI1:barkI2)...
                                         + earDFAdjust(barkI1:barkI2);
            case "noOuter"
                excitation = lvlWeightCB;
        end
    else
        excitation = lvlWeightCB;
    end
    
    % transformation to loudness
    maskLTQ = excitation > levelThresQ(barkI1:barkI2);  % mask used later to zero below-threshold excitation
    excitCBands = excitation - critBWAdjust(barkI1:barkI2);
    levelThresQRep = repmat(levelThresQ(barkI1:barkI2), size(excitation, 1), 1);

    switch adjustLoud
        case 'iso5321'
            loudCore = calN.*((1 - threshFact...
                               + threshFact.*10.^((excitCBands...
                                                   - levelThresQRep)/10)).^0.25...
                              - 1).*0.0635.*(10.^(0.1*levelThresQRep)).^0.25;

        case {'iso5323', 'ecma4182'}
        
            if strcmp(adjustLoud, 'ecma4182')
                tFormParamsCBands = tFormParams(:, barkI1:barkI2);

            else
                tFormParamsCBands = tFormParams;
            end

            tFormP1 = tFormParamsCBands(1, :);
            tFormP2 = tFormParamsCBands(2, :);
            tFormP3 = tFormParamsCBands(3, :);
            tFormP4 = tFormParamsCBands(4, :);
            tFormP5 = tFormParamsCBands(5, :);
            tFormP6 = tFormParamsCBands(6, :);
            tFormP7 = tFormParamsCBands(7, :);
            
            tFormExpt = tFormP2.*exp(-tFormP3.*excitCBands) + tFormP4...
                        + tFormP5.*tanh((excitCBands - tFormP6)./tFormP7);
    
            loudCore = calN.*((1 - threshFact...
                               + threshFact.*10.^((excitCBands...
                                                   - levelThresQRep)/10)).^tFormExpt...
                              - 1).*tFormP1.*(10.^(0.1*levelThresQRep)).^tFormExpt;

    end

    loudCore(~maskLTQ) = 0;  % zero loudness for below-threshold excitation
    loudCore(loudCore < 0) = 0;  % zero negative loudness

    % correction to core loudness in lowest critical band (Annex
    % A.4 f_corr_loudness) "for the consideration of the run of threshold
    % in quiet within this critical band"
    if CB1
        corrCL = (0.4 + 0.32.*loudCore(:, 1).^0.2);
        corrCL(corrCL >= 1) = 1;
        loudCore(:, 1) = loudCore(:, 1).*corrCL;
    end

    % Specific loudness slopes adjustment
    % add dummy band for uppermost band
    if barkNOut(end) >= 23.6
        loudCore1 = [loudCore, zeros(size(loudCore, 1), 1)];
    else
        loudCore1 = loudCore;
    end
    
    %%%% SQAT(/AARAE) code adapted
    
    loudTDepChan = zeros(nTimeSteps, 1);
    specLoudChan = zeros(nTimeSteps, barkAxisN);
    
    for l = 1:nTimeSteps
    
        N = 0;
        z1 = 0;  % critical band rate starts at 0
        n1 = 0;  % loudness level starts at 0
        iz = 1;
        z = 0.1;
        j = 18;
    
        for i = 1:barkCount + 1  % specific loudness
    
            % Determines where to start on the slope
            ig = i - 1;
    
            % steepness of upper slope (specNSlope) for bands above 8th one are identical
            if ig > 8
                ig = 8;
            end
    
            while z1 < barkN(i)
    
                if n1 <= loudCore1(l, i)     % loudCore1 is the main loudness level
                    % contribution of unmasked main loudness to total loudness
                    % and calculation of values
                    if n1 < loudCore1(l, i)
                        j = 1;
    
                        while (specNSlopeBound(j) > loudCore1(l, i)) && (j < 18) % the value of j is used below to build a slope
                            j = j + 1; % j becomes the index at which Nm(i)                        % to the range of specific loudness
                        end
                    end
    
                    z2 = barkN(i);
                    n2 = loudCore1(l, i);
                    N = N + n2*(z2 - z1);
                    k = z;                     % initialisation of k
    
                    while (k <= z2)
                        specLoudChan(l, iz) = n2;
                        iz = iz + 1;
                        k = k + (1/10);
                    end
    
                    z = k;
    
                else %if n1 > loudCore(i)
                    % decision wether the critical band in question is completely
                    % or partly masked by accessory loudness
    
                    n2 = specNSlopeBound(j);
    
                    if n2 < loudCore1(l, i)
                        n2 = loudCore1(l, i);
                    end
    
                    dz = (n1 - n2) / specNSlope(j, ig);
                    z2 = z1 + dz;
    
                    if z2 > barkN(i)
                        z2 = barkN(i);
                        dz = z2 - z1;
                        n2 = n1 - dz*specNSlope(j, ig);
                    end
    
                    N = N + dz*(n1 + n2)/2;
                    k = z;                     % initialisation of k
    
                    while (k <= z2)
                        specLoudChan(l, iz) = n1 - (k - z1)*specNSlope(j, ig);
                        iz = iz + 1;
                        k = k + (1/10);
                    end
    
                    z = k;
    
                end
    
                if (n2 <= specNSlopeBound(j)) && (j < 18)
                    j = j + 1;
                end
    
                if (n2 <= specNSlopeBound(j)) && (j >= 18)
                    j = 18;
                end
    
                z1 = z2;     % n1 and z1 for next loop
                n1 = n2;
    
            end
        end
    
        if N < 0
            N = 0;
        end
    
        if N <= 16
            N = (N*1000 + 0.5)/1000;
        else
            N = (N*100 + 0.5)/100;
        end
    
        loudTDepChan(l) = N; % total loudness at current timeframe l
    end
    
    %%%% End of SQAT(/AARAE) code adapted
    specLoud(:, :, chan) = specLoudChan;
    loudTDep(:, chan) = loudTDepChan;
end  % end of for loop over channels

% power-averaged overall loudness
loudPowAvg = power(mean(loudTDep.^(1/log10(2)), 1), log10(2));

% 95th percentile overall loudness
loud5pcEx = prctile(loudTDep, 95, 1);

% time-averaged specific loudness as a function of Bark number
specLoudAvg = squeeze(mean(specLoud, 1));

% loudness level
loudLevel = 40*(loudTDep + 5e-5).^0.35;
loudLevel(loudTDep >= 1) = 40 + 10*log2(loudTDep(loudTDep >= 1));

if loudPowAvg < 1
    loudLvlPowAvg = 40*(loudPowAvg + 5e-5).^0.35;
    loudLvl5pcEx = 40*(loud5pcEx + 5e-5).^0.35;
else
    loudLvlPowAvg = 40 + 10*log2(loudPowAvg);
    loudLvl5pcEx = 40 + 10*log2(loud5pcEx);
end

% assign outputs
loudness.loudTDep = loudTDep;
loudness.loudPowAvg = loudPowAvg;
loudness.loud5pcEx = loud5pcEx;
loudness.loudLevel = loudLevel;
loudness.loudLvlPowAvg = loudLvlPowAvg;
loudness.loudLvlPowAvg = loudLvl5pcEx;
loudness.specLoudAvg = specLoudAvg;
loudness.specLoud = specLoud;
loudness.barkAxis = barkAxis;
loudness.soundField = soundField;
loudness.adjustLoud = adjustLoud;
loudness.freqInMid = freqInMid;
loudness.freqInNom = freqInNom;

end  % end of acousticQuasiLoudZwicker function