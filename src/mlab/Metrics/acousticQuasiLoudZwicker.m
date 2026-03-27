function loudness = acousticQuasiLoudZwicker(spectrL, fLim, timeStep, axisF, soundField, adjustLoud, outPlot)
% loudness = acousticQuasiLoudZwicker(spectrL, fLim, timeStep, axisF,
%                                     soundField, adjustLoud, outPlot)
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
%   Contiguous input sound level 1/3-octave Leq spectrum or spectra for
%   processing. Must be orientated with time and frequency bands on
%   first two axes, and channels on third axis.
%
% fLim : vector (default: [25, 12600])
%   Frequency limits for the input spectra (these are
%   automatically matched to nearest exact band centre-frequencies).
%
% timeStep : number (default: 0.1)
%   Time step value (seconds) used for time-dependent inputs.
%
% axisF : integer optional (1 or 2, default: 2)
%   Fequency band axis for series of spectra.
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
% outPlot : Boolean true/false (default: false)
%   Flag indicating whether to generate a figure from the output.
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
% soundField : string
%   Identifies the soundfield type applied (the input argument soundField)
%
% timeOut : vector
%   Time (seconds) corresponding with time-dependent outputs.
%
% adjustLoud : string
%   Indicates which loudness model standard was adjusted for (=adjustLoud
%   input).
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
% Paulus, E & Zwicker, E, 1972. 
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
% Date last modified: 27/03/2026
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
        fLim (1, 2) double {mustBeInRange(fLim, 25, 12600)} = [25, 12600]
        timeStep (1, 1) double {mustBePositive} = 0.1
        axisF (1, 1) {mustBeInteger, mustBeInRange(axisF, 1, 2)} = 2
        soundField (1, :) string {mustBeMember(soundField,...
                                               {'freeFrontal',...
                                                'diffuse',...
                                                'noOuter'})} = 'freeFrontal'
        adjustLoud (1, :) string {mustBeMember(adjustLoud,...
                                               {'iso5321',...
                                                'iso5323',...
                                                'ecma4182'})} = 'iso5321'
        outPlot (1, 1) {mustBeNumericOrLogical} = false
    end

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("src", "mlab")))

%% Input checks
if ndims(spectrL) == 3
    inChans = size(spectrL, 3);
elseif ismatrix(spectrL)
    inChans = 1;
else
    error("Error: Input spectrL must not have more than 3 dimensions.")
end

if axisF == 1
    if inChans == 1
        spectrL = spectrL.';
    else
        spectrL = permute(spectrL, [2, 1, 3]);
    end
end

if size(spectrL, 1)*timeStep < 0.5
    error("Error: the input signal duration is too short (it should be at least 0.5 s).")
end

% fractional octave frequencies
[freqInMid, freqLower, freqUpper, freqInNom] = noctf(fLim, 3);
[fmAll, ~, fAllUpper, ~] = noctf([25, 12600], 3);

if length(freqInMid) ~= size(spectrL, 2)
    error("Error: The frequency band range of the input spectra must correspond with the input band limits 'fLim'.")
end

%% Define constants

% frequency band limiting indices
fmAllI1 = find(min(freqInMid) == fmAll);
fmAllI2 = find(max(freqInMid) == fmAll);
fmShift1 = fmAllI1 - 1;

% number of time steps in spectral series
nTimeSteps = size(spectrL, 1);

% ISO 532-1 Table A.3: low-frequency (25-250 Hz) weights
lowFWeights = [-32, -24, -16, -10, -5, 0, -7, -3, 0, -2, 0;
               -29, -22, -15, -10, -4, 0, -7, -2, 0, -2, 0;
               -27, -19, -14,  -9, -4, 0, -6, -2, 0, -2, 0;
               -25, -17, -12,  -9, -3, 0, -5, -2, 0, -2, 0;
               -23, -16, -11,  -7, -3, 0, -4, -1, 0, -1, 0;
               -20, -14, -10,  -6, -3, 0, -4, -1, 0, -1, 0;
               -18, -12,  -9,  -6, -2, 0, -3, -1, 0, -1, 0;
               -15, -10,  -8,  -4, -2, 0, -3, -1, 0, -1, 0];

lowFWeightsAdj = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

[lfwRows, lfwCols] = size(lowFWeights);

switch adjustLoud
    case 'iso5323'     % 25 31.5 40  50  63  80 100 125 160 200 250
        lowFWeightsAdj = [30, 27, 21, 22, 20, 18, 5, 3, 2, -2, -3];
                       % 25 31.5 40  50  63  80 100 125 160 200 250
        lowFweightsAdj2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                          -1, -2, 0, -1, -1, -1, -1, -1, 0, 0, 0;
                          -2, -5, -2, -3, -4, -3, -1, -1, 0, 0, 0;
                          -2, -5, -2, -4, -5, -3, -1, -1, 0, 0, 0;
                          -7, -8, -3, -7, -8, -7, -4, -4, -1, 0, 0;
                          -7, -10, -6, -9, -9, -10, -4, -4, -1, 0, 0;
                          -9, -12, -7, -10, -11, -13, -2, -2, -2, 0, 0;
                          -9, -12, -8, -12, -12, -13, -2, -2, -2, 0, 0];

        lowFWeights = lowFWeights + lowFWeightsAdj + lowFweightsAdj2;

    case 'ecma4182'    % 25 31.5 40  50  63  80 100 125 160 200 250
        lowFWeightsAdj = [42, 36, 27, 24, 22, 19, 11, 9, 7, 0, -2];
                        % 25 31.5 40  50  63  80 100 125 160 200 250
        lowFweightsAdj2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                          -3, -4, 0, -1, -1, -1, -1, -1, 0, 0, 0;
                          -5, -6, -1, -2, -3, -2, -2, -2, 0, 0, 0;
                          -9, -10, -3, -4, -5, -3, -3, -2, 0, 0, 0;
                          -13, -13, -6, -7, -8, -6, -5, -5, -1, 0, 0;
                          -17, -16, -10, -10, -10, -9, -5, -5, -1, 0, 0;
                          -20, -19, -12, -11, -12, -11, -5, -3, -2, 0, 0;
                          -20, -19, -13, -13, -13, -11, -5, -3, -2, 0, 0];

        lowFWeights = lowFWeights + lowFWeightsAdj + lowFweightsAdj2;

end

% insert zeros to match all bands and tranpose to ease later calcs
lowFWeights = [lowFWeights, zeros([lfwRows, length(fmAll) - lfwCols])].';

lowFWeightsAdj = [lowFWeightsAdj, zeros([1, length(fmAll) - lfwCols])];

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

% adjustment for original Zwicker third-octave definitions to 
% critical band limits from Paulus et al, 1972, Table II
% fgZ = [0, 90, 180, 280, 355, 450, 560, 710, 900, 1120, 1400,...
%       1800 2240 2800 3550 4500 5600 7100 9000 11200 14000];

% dfgZ = diff(fgZ);

% mid-frequencies for frequency bands corresponding with critical
% bandwidths used to determine ISO 532-1 Table A.7, from Paulus et al,
% 1972, Table II
fmZ = [63, 125, 224, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000,...
       2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500];

% Zwicker approximation for critical bandwidth
cbwZ = 25 + 75*(1 + 1.4*(fmZ/1000).^2).^0.69;

% reproduction of data points in Paulus et al, 1972, Figure 5
% (through which a curve was fit to smooth the bandwidth adjustments)
% dLTG = 10*log10(dfgZ./cbwZ);

% implied bandwidths obtained from the smoothed adjustment curve
fbwZ = cbwZ.*10.^(critBWAdjust/10);

% third-octave bandwidths corresponding with critical bands
fbw = [fAllUpper(6), fAllUpper(9) - fAllUpper(6), fAllUpper(11) - fAllUpper(9), diff(fAllUpper(11:end))];

% level adjustment corresponding with difference between Zwicker
% third-octave bandwidths and present bandwidths
dLfbw = 10*log10(fbw./fbwZ);

switch adjustLoud
    case 'iso5321'
        % ISO 532-1 Table A.6: critical band level at threshold in quiet
        levelThresQ = [30, 18, 12, 8, 7, 6, 5, 4,...
                       3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3];

    case 'iso5323'
        % ISO 523-3 Table 2 (extended to 25 Hz as per ISO 532-3 code)
        levelThresQTOB = [28.18, 28.18, 28.18, 28.18, 23.9, 19.2, 15.68,...
                          12.67, 10.09, 8.08, 6.3, 5.3, 4.5, 3.63, 3.63,...
                          3.63, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63,...
                          3.63, 3.63, 3.63, 3.63, 3.63];

        levelThresQCB = [10*log10(sum(10.^(levelThresQTOB(1:6)/10))),...
                         10*log10(sum(10.^(levelThresQTOB(7:9)/10))),...
                         10*log10(sum(10.^(levelThresQTOB(10:11)/10))),...
                         levelThresQTOB(12:end)];

        levelThresQ = levelThresQCB - critBWAdjust + dLfbw;

    case 'ecma4182'
        % ISO 523-3 Table 2 (extended to 25 Hz with extrapolation)    
        levelThresQTOB = [45, 39, 33, 28.18, 23.9, 19.2, 15.68,...
                          12.67, 10.09, 8.08, 6.3, 5.3, 4.5, 3.63, 3.63,...
                          3.63, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63,...
                          3.63, 3.63, 3.63, 3.63, 3.63];

        levelThresQCB = [10*log10(sum(10.^(levelThresQTOB(1:6)/10))),...
                         10*log10(sum(10.^(levelThresQTOB(7:9)/10))),...
                         10*log10(sum(10.^(levelThresQTOB(10:11)/10))),...
                         levelThresQTOB(12:end)];

        levelThresQ = levelThresQCB - critBWAdjust + dLfbw;

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

% calculate number of critical bands for output
barkNMapOut = barkNMap(fmAllI1:fmAllI2);
barkNOut = unique(barkNMapOut);
barkCount = length(barkNOut);
barkI1 = find(min(barkNOut) == barkN);
barkI2 = find(max(barkNOut) == barkN);

if barkI1 == 1
    barkMin = 0.1;
else
    barkMin = barkN(barkI1 - 1) + 0.1;
end
barkMax = max(round(barkN(barkI2 + 1), 1));
barkAxis = barkMin:0.1:barkMax;
barkAxisN = length(barkAxis);
barkNRange = barkN(barkI1:barkI2 + 1);  % redefines bark numbers according to critical band range selected

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
% 40 dB in free-field, when using ISO 532-1 1/3-octave filters)
switch adjustLoud
    case 'iso5321'
        calN = 1.214117275158246;%/1.003122045727327/0.999417805029527/1.000108620697671;
    case 'iso5323'
        calN = 1.286531792122494;%/1.144265905393895/1.009638100259755/0.998534527407563/1.000224091277673;
    case 'ecma4182'
        calN = 1.129367701604459;%/1.109032916997183/0.992236701701909/1.000077784186351;
end

% loudness transform threshold factor
threshFact = 0.25;

% loudness transform parameters (per critical band)
switch adjustLoud
    case 'iso5323'
        tFormParams = [0.005568892, 0.034082322, 0.05276281, 0.047589615,...
                       0.055233763, 0.060735963, 0.066341805, 0.059079492,...
                       0.051050047, 0.0506713, 0.053853145, 0.054256423,...
                       0.055192073, 0.060253342, 0.050565467, 0.063925214,...
                       0.063033356, 0.084145715, 0.090950214, 0.10455592;
                       -0.25, -0.25, -0.250000138, -0.250020838,...
                       -0.301393576, -0.35179706, -0.25508033, -0.340899753,...
                       -0.25, -0.324658385, -0.276203846, -0.257634047,...
                       -0.15, -0.1425, -0.1357, -0.24,...
                       -0.374591987, -0.250148994, -0.28650851, -0.244853157;
                       0.5, 0.750000843, 0.748949038, 0.760869545,...
                       0.092128909, 0.098124703, 0.088494106, 0.1,...
                       0.100003056, 0.112069958, 0.107052359, 0.076757664,...
                       0.089107633, 0.105400936, 0.106375595, 0.081944186,...
                       0.096786354, 0.055322228, 0.062490766, 0.066315428;
                       0.343812576, 0.304689875, 0.285528292, 0.290264531,...
                       0.283415709, 0.276547317, 0.273826644, 0.278848845,...
                       0.29394038, 0.284097732, 0.282592126, 0.279969074,...
                       0.2779968, 0.272629845, 0.279427994, 0.312993714,...
                       0.167764298, 0.24570819, 0.223082005, 0.026073762;
                       -0.065695308, -0.014665965, -0.0157655, -0.012151657,...
                       -0.013923259, -0.014076095, -0.013380292, -0.011307706,...
                       -0.021332485, -0.00843291, -0.008336012, -0.005874745,...
                       -0.003568029, -0.004715027, -0.000400398, -0.046745988,...
                       9.41E-02, -6.17E-04, 1.36E-02, 0.190602182;
                       -5.873556067, 101.3580165, 100.4504578, 100.6015338,...
                       100.0243422, 100.452738, 100.156096, 100.3243,...
                       119.9975644, 104.1265168, 106.2842948, 102.2741459,...
                       99.33557688, 99.00725743, 2.23352E-05, 121.4117389,...
                       161.4909332, 191.5686578, 214.7061278, 113.5427639;
                       -174.8631934, -24.60374535, -23.40908683, -20.08079917,...
                       -21.69198323, -22.20188116, -21.86075832, -20.83095038,...
                       -43.10151071, -20.45443513, -22.54450572, -17.80208162,...
                       -13.21897145, -23.86138725, -59.06859625, -3.386149498,...
                       -6.117613803, -0.002056314, -81.70886929, -0.682841422;
                       0, 0.025, 0.02, 0,...
                       0.005, 0.00278, 0.0025, 0.0018,...
                       0.0015, 0.002, 0.00115, 0,...
                       0.00185, 0.00001, 0.00125, 0,...
                       0.002, 0, 0.005, 0.026];

    case 'ecma4182'
        tFormParams = [0.193940725, 0.22473418, 0.190592035, 0.112106116,...
                       0.101291184, 0.109015913, 0.109374371, 0.118245569,...
                       0.109999749, 0.119372397, 0.112684242, 0.124611825,...
                       0.113166473, 0.106665524, 0.102058517, 0.115422121,...
                       0.127086523, 0.144206986, 0.158044149, 0.160612695;
                       81.29436687, 0.312958826, 0.5, -0.25,...
                       -0.249999988, -7.36073E-05, -0.249999931, -0.250008575,...
                       1.000021114, 12.5223323, 0.231344739, -0.099194902,...
                       -0.00105586, 0, 0, 0,...
                       4.562278532, 10.46361192, 11.56033451, -0.000385932;
                       0.131730766, 0.058585588, 0.2, 1.040040712,...
                       0.842940403, 0.000117904, 1.352234593, 0.625075741,...
                       0.600000098, 1.010756852, 0.458497745, 0.712716374,...
                       0.800647094, 0.538775212, 0.890175983, 0.750378202,...
                       0.827418154, 0.813464016, 0.76962153, 0.9446303;
                       0.084916156, 0.025659201, 0.167546059, 0.19858002,...
                       0.205674409, 0.207375617, 0.204649396, 0.206616452,...
                       0.209885744, 0.203529183, 0.205523548, 0.229430347,...
                       0.196796311, 0.198900145, 0.195802129, 0.19230645,...
                       0.191230317, 0.187937737, 0.182029422, 0.172055865;
                      -0.184958038, -0.22197563, -0.049786187, -0.029460114,...
                      -0.018156714, -0.025388855, -0.021304891, -0.029243358,...
                      -0.022235612, -0.021318902, -0.021067779, -0.052628765,...
                      -0.010881264, -0.009874284, -0.004979137, -0.010540424,...
                      -0.017459012, -0.025178506, -0.028558922, -0.02386538;
                       2.71288E-15, -0.117532982, 92.55691885, 99.73825631,...
                       104.5831917, 115.9862055, 109.8230144, 120.7912931,...
                       119.1148516, 116.2533806, 115.2210912, 143.1272494,...
                       106.9099775, 104.2084333, 99.82472834, 106.0163686,...
                       110.5616377, 111.9317432, 109.9240041, 104.319429;
                      -200, -138.6192611, -62.71809962, -54.62700651,...
                      -37.13555077, -49.76827964, -41.55740267, -49.35390318,...
                      -35.06377794, -33.02100315, -33.24526782, -43.1501127,...
                      -17.65301799, -16.05886478, -2.360707924, -13.61587154,...
                      -20.04941176, -25.2808917, -25.99735934, -20.34042365;
                      0.00071758, 0, 0.018, 0,...
                      0.007052042, 0.011937378, 0.014484559, 0.0125,...
                      0.018000276, 0.004994058, 0.015, 0.009924449,...
                      0.014, 0.013, 0.0125, 0,...
                      0, -0.015, -0.015, -0.015];

end  % end of adjustLoud switch

% Core loudness parameters for scaling adjustments
switch adjustLoud
    case 'iso5321'
        loudCoreScalLim = [0.175206559279323, 0.0998721141962912,...
                           0.0562259655388907, 0.0503291994566649,...
                           0.0323787742019847, 0.0214724333197948,...
                           0.0145663061860690, 0.0106963917147328,...
                           0.0120397538551626, 0.0156867555722016,...
                           0.0145351727954897, 0.00792706152583279,...
                           0.00530757169293728, 0.00716388421756555,...
                           0.0108774288315528, 0.0188497593204034,...
                           0.0438528596384881, 0.0692933017387570,...
                           0.0714963393585361, 0.0211255209023633;
                           40.8915609307780, 33.1919723833486,...
                           31.1766462873821, 30.5157031895724,...
                           29.7440782185720, 28.8317267592538,...
                           28.0752521555016, 27.6155419727084,...
                           28.0406253119820, 32.9098829123404,...
                           37.0788723352572, 33.4449633499979,...
                           31.7412262374187, 34.5193582730931,...
                           38.0930314504745, 42.5129104626755,...
                           49.2861715403145, 50.0674042739489,...
                           40.2310227614464, 16.7046470437925];

        polyCoeffs = [-0.0177267146309573, -0.0181524729759670, -0.0151004824102610,...
                      -0.0204504698683761, -0.0145797277413292, -0.00970899369407536,...
                      -0.00658465286305260, -0.00658951468042123, -0.00658951468042123,...
                      -0.00795199012848177, -0.00681059843461716, -0.00622253857462985,...
                      -0.00553819307861815, -0.00614570719493355, -0.00777665965838361,...
                      -0.0120768951375648, -0.0157452223046718, -0.0150301298430959,...
                      -0.0121384793351893, -0.00499525056758499;
                      0.103795799922857, 0.103995861300837, 0.0958173900261008,...
                      0.115770875297616, 0.0967247295950414, 0.0829120073515128,...
                      0.0703988379961275, 0.0676721783299202, 0.0676721783299202,...
                      0.0658343307032811, 0.0602722360455543, 0.0592923883846582,...
                      0.0606031653241484, 0.0599424497156373, 0.0654941735945591,...
                      0.0771629715172939, 0.0840167915075612, 0.0731189240417204,...
                      0.0515867874544800, 0.0134739166630018;
                      1.13910918105779, 1.12017536062933, 1.10744166940714,...
                      1.10776080800244, 1.10441452052549, 1.10100662082541,...
                      1.10797231061891, 1.11850630410739, 1.11850630410739,...
                      1.12871113619910, 1.13198477679486, 1.12855885886580,...
                      1.12075339775294, 1.12347320049110, 1.12252020177000,...
                      1.12934585209696, 1.15607058604635, 1.18159147322955,...
                      1.20057521209692, 1.22434987865123];

        xShift = 0;  % null shift for ISO 532-1 log polynomial fit

    case 'iso5323'
        loudCoreScalLim = [0.0686232877960075, 0.0455369521558487,...
                           0.0304510825882214, 0.0239249691738821,...
                           0.0189775994866060, 0.0138938815719889,...
                           0.0136962475253381, 0.00876840611351073,...
                           0.00915593829481202, 0.00969672193834636,...
                           0.0119378796574062, 0.00843887026528913,...
                           0.0155180247622181, 0.0152055682075589,...
                           0.0166645927344482, 0.0125882987565647,...
                           0.0122058046979527, 0.0104817460307530,...
                           0.0194669677922798, 0.0463091733969044;
                           45.0768226453848, 37.5149169006044,...
                           37.2111048886061, 46.0782975233605,...
                           52.5681668963248, 52.7595117394585,...
                           57.8711048286698, 59.1305325223868,...
                           58.6333613989231, 63.6086968236880,...
                           79.9530502749036, 77.7106476101222,...
                           76.2163483675857, 76.2744728746143,...
                           75.4673426532020, 71.8137199984029,...
                           53.6698765772943, 39.4152272851462,...
                           29.1455836371888, 14.7998939008256];

        polyCoeffs = [0.0227931208305138, 0.0285738404030923, 0.00541711466401674,...
                      0.00990169222264163, 0.0145130279348129, 0.00710242135098269,...
                      0.0130651869878324, 0.00352262970626301, 0.00352262970626301,...
                      0.00664536668191815, 0.0183814729003890, 0.00185106230660431,...
                      0.00258041914559566, -0.00550840826085959, -0.0212048809446684,...
                      -0.0324808223367526, -0.0556344698253639, -0.0698773692357039,...
                      -0.0482012781527130, 0.0513332689503584;
                      0.351298099668721, 0.138852122388597, 0.258883772621513,...
                      0.323743726246527, 0.330715914123338, 0.372837142871287,...
                      0.342732073161703, 0.391807933953650, 0.391807933953650,...
                      0.371231349201823, 0.311908987065506, 0.406554582433415,...
                      0.360537801521611, 0.376299108872329, 0.456838271845548,...
                      0.506828905862457, 0.573760846253329, 0.598913305595168,...
                      0.411826268687735, -0.144304637397255;
                      0.655644460064169, 0.786030981577344, 0.698158646389966,...
                      0.563992449646277, 0.640691152330000, 0.667158665039848,...
                      0.755372087380492, 0.703650705897034, 0.703650705897034,...
                      0.771323765595285, 0.873351461265575, 0.695525773651131,...
                      0.774716056400417, 0.838667209193858, 0.700960581335911,...
                      0.637932135941241, 0.775290679808694, 0.694493801393624,...
                      0.605110968916654, 0.896748225671540];

        xShift = 1;  % shift for ISO 532-3 log polynomial fit

    case 'ecma4182'
        loudCoreScalLim = [0.134120067887792, 0.0709563870814197,...
                           0.0339494157810020, 0.0273280219437117,...
                           0.0320310575778420, 0.0383482371619014,...
                           0.0393956633238017, 0.0363203408289930,...
                           0.0375295594340747, 0.0181995647038982,...
                           0.0295908274847412, 0.0231527609564899,...
                           0.0257160154785820, 0.0247027241319666,...
                           0.0300141631463563, 0.0224273147233404,...
                           0.0267222571126299, 0.0154566667107326,...
                           0.0233693354099823, 0.0117780726621573;
                           13.2601893858513, 13.8039143413788,...
                           10.8919144283017, 12.1202503699032,...
                           13.3445900992356, 14.0952338589571,...
                           15.3681216900662, 15.5856698208251,...
                           15.0200397781887, 14.6336278043039,...
                           17.1799810026469, 16.3139494752907,...
                           15.7392465853213, 17.0169367190822,...
                           17.0577577128619, 19.1431915170042,...
                           18.0108441557629, 15.8337724033614,...
                           13.4493814168596, 6.98750705535467];

        polyCoeffs = [0.0281813130454864, 0.0398579909047190, -0.0311067706969822,...
                      -0.0579178421063959, -0.0430899545094520, -0.0307751333344124,...
                      -0.0338038378518385, -0.0487986324881162, -0.0487986324881162,...
                      -0.109907011369036, -0.0665332895465800, -0.0756531933618948,...
                      -0.0647363721767976, -0.0656845805198346, -0.0889903776037637,...
                      -0.106801885650847, -0.118328603793767, -0.166732106691452,...
                      -0.145663974957413, -0.114753616319905;
                      0.0999134759157014, 0.163879680024510, 0.413988340468798,...
                      0.543590862892810, 0.494938536117935, 0.444855579491137,...
                      0.455273354431365, 0.512417725042014, 0.512417725042014,...
                      0.791788485868339, 0.579329031009410, 0.642541596561724,...
                      0.589890603633805, 0.578333087885605, 0.661602100419081,...
                      0.748375468354315, 0.781146248978933, 0.940754983264441,...
                      0.771489073512790, 0.547505859332814;
                      1.07620218014001, 0.934721289410464, 0.632173325143473,...
                      0.405288470772238, 0.534426827140865, 0.651211779719732,...
                      0.734135624114756, 0.711714336601515, 0.711714336601515,...
                      0.670056548726840, 0.918155480292106, 0.725083400560378,...
                      0.727573118591267, 0.770367505406647, 0.679598092224904,...
                      0.593865643218321, 0.744668788744825, 0.575812678137361,...
                      0.435968156429125, 0.375756434108614];

        xShift = 1;  % shift for ECMA-418-2 log polynomial fit
end


%% Signal processing

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
        if any(lowFWeights(kk, :)) %&& strcmp(adjustLoud, 'iso5321')
    
            % repeat band levels and apply all weights
            spectrLRepWt = repmat(spectrL(:, kk, chan), 1, rangeInts) + lowFWeights(kk, :);
    
            % check which weighted levels meet corresponding range criterion
            meetCriterion = (spectrLRepWt - lowFWeightsAdj(kk)) <= lowFLRange;
    
            % get column indices for lowest criterion match
            [~, lowCriterionMatch] = max(meetCriterion == 1, [], 2);
    
            % repeat weights matrix for each time step
            lowFWeightsRep = repmat(lowFWeights(kk, :), nTimeSteps, 1);
    
            % weight levels according to level criteria
            levelWeighted(:, kk) = spectrL(:, kk, chan)...
                                   + lowFWeightsRep(sub2ind(size(lowFWeightsRep),...
                                   (1:nTimeSteps).', lowCriterionMatch));
        end  % end of if branch for non-zero weighting
    end  % end of for loop through bands
    
    % apply if using ECMA or ISO 532-3 outer-middle ear filter functions
    if ~strcmp(adjustLoud, 'iso5321')
        switch soundField
            case 'freeFrontal'
                levelWeighted = levelWeighted + outMidFree(fmAllI1:fmAllI2);
            case 'diffuse'
                levelWeighted = levelWeighted + outMidDiffuse(fmAllI1:fmAllI2);
            case 'noOuter'
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
    if any(barkNMapOut == barkN(2)) && fmAllI2 >= 9
        
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
    if any(barkNMapOut == barkN(3)) && fmAllI2 >= 11
    
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
    % apply ISO 532-1 outer ear filter
    if strcmp(adjustLoud, 'iso5321')
        switch soundField
            case 'freeFrontal'
                excitation = lvlWeightCB - earTmission(barkI1:barkI2);
            case 'diffuse'
                excitation = lvlWeightCB - earTmission(barkI1:barkI2)...
                                         + earDFAdjust(barkI1:barkI2);
            case 'noOuter'
                excitation = lvlWeightCB;
        end
    else
        % Third-octave bandwidth adjustment
        excitation = lvlWeightCB + dLfbw(barkI1:barkI2);
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

            tFormParamsCBands = tFormParams(:, barkI1:barkI2);

            tFormP1 = tFormParamsCBands(1, :);
            tFormP2 = tFormParamsCBands(2, :);
            tFormP3 = tFormParamsCBands(3, :);
            tFormP4 = tFormParamsCBands(4, :);
            tFormP5 = tFormParamsCBands(5, :);
            tFormP6 = tFormParamsCBands(6, :);
            tFormP7 = tFormParamsCBands(7, :);
            tFormP8 = tFormParamsCBands(8, :);

            tFormExpt = tFormP2.*exp(-tFormP3.*excitCBands) + tFormP4...
                        + tFormP5.*tanh((excitCBands - tFormP6)./tFormP7);

            loudCore = calN.*(((1 - threshFact...
                               + threshFact.*10.^((excitCBands...
                                                   - levelThresQRep)/10)).^tFormExpt...
                              - 1).*tFormP1.*(10.^(0.1*levelThresQRep)).^tFormExpt...
                              + tFormP8);

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

    % setup for scaling adjustments
    loudCoreScalVals = loudCore;
    loudCoreScalMin = repmat(loudCoreScalLim(1, :), [size(loudCore, 1), 1]);
    maskMin = loudCore < loudCoreScalMin;
    loudCoreScalVals(maskMin) = loudCoreScalMin(maskMin);
    loudCoreScalMax = repmat(loudCoreScalLim(2, :), [size(loudCore, 1), 1]);
    maskMax = loudCore > loudCoreScalMax;
    loudCoreScalVals(maskMax) = loudCoreScalMax(maskMax);
    loudCoreScalTform = log2(loudCoreScalVals + xShift);
    
    % loop through critical bands and estimating scaling factors
    loudCoreScale = ones(size(loudCoreScalTform));
    for iBand = 1:barkCount
        loudCoreScale(:, iBand) = polyval(polyCoeffs(:, iBand), loudCoreScalTform(:, iBand));
    end

    % these operations reduce the effect of the scaling at low sone
    % values, and force it to unity at 1 sone (~0.4 sone core loudness).
    maskCompress = round(loudCore, 1) <= 0.4;
    loudCoreScale(maskCompress) = (loudCoreScale(maskCompress)).^(1/4);

    maskSone = round(loudCore, 1) == 0.4;
    loudCoreScale(maskSone) = 1;

    % apply scaling
    loudCoreScaled = loudCore./loudCoreScale;

    % Specific loudness slopes adjustment
    % add dummy band for uppermost band
    loudCore1 = [loudCoreScaled, zeros(size(loudCoreScaled, 1), 1)];

    %% SQAT(/AARAE) code adapted
    
    loudTDepChan = zeros(nTimeSteps, 1);
    specLoudChan = zeros(nTimeSteps, barkAxisN);
    
    for lStep = 1:nTimeSteps
    
        N = 0;
        n1 = 0;  % loudness starts at 0
        if CB1
            z0 = 0;  % critical band rate starts at 0
        else
            z0 = min(barkNRange);
        end
        z1 = z0 + 0.1;
        iz = 1;  % critical band index
        jNidx = 18;  % index of loudness for slope calculation (Table A.9)
    
        for iBand = 1:barkCount + 1  % specific loudness
    
            % Determines where to start on the slope
            ig = iBand - 1 + (barkI1 - 1);  % the addition accounts for a start band above band 1
    
            % steepness of upper slope (specNSlope) for bands above 8th one are identical
            if ig > 8
                ig = 8;
            end
    
            while z0 < barkNRange(iBand)
    
                if n1 <= loudCore1(lStep, iBand)     % loudCore1 is the core loudness
                    % contribution of unmasked main loudness to total loudness
                    % and calculation of values
                    if n1 < loudCore1(lStep, iBand)
                        jNidx = 1;
    
                        % the value of jNidx is used below to build a slope to the range of specific loudness
                        while (specNSlopeBound(jNidx) > loudCore1(lStep, iBand)) && (jNidx < 18)
                            jNidx = jNidx + 1; % jNidx becomes the index at which Nm(iBand)
                        end
                    end
    
                    z2 = barkNRange(iBand);
                    n2 = loudCore1(lStep, iBand);
                    N = N + n2*(z2 - z0);
                    zk = z1;                     % initialisation of zk
    
                    while (zk <= z2)
                        specLoudChan(lStep, iz) = n2;
                        iz = iz + 1;
                        zk = zk + 0.1;
                    end
    
                    z1 = zk;
    
                else  % if slope loudness > core loudness
                    
                    % decision whether the critical band in question is completely
                    % or partly masked by accessory loudness
    
                    n2 = specNSlopeBound(jNidx);
    
                    if n2 < loudCore1(lStep, iBand)
                        n2 = loudCore1(lStep, iBand);
                    end
    
                    dz = (n1 - n2) / specNSlope(jNidx, ig);
                    z2 = z0 + dz;
    
                    if z2 > barkNRange(iBand)
                        z2 = barkNRange(iBand);
                        dz = z2 - z0;
                        n2 = n1 - dz*specNSlope(jNidx, ig);
                    end
    
                    N = N + dz*(n1 + n2)/2;
                    zk = z1;                     % initialisation of zk
    
                    while (zk <= z2)
                        specLoudChan(lStep, iz) = n1 - (zk - z0)*specNSlope(jNidx, ig);
                        iz = iz + 1;
                        zk = zk + (1/10);
                    end
    
                    z1 = zk;
    
                end
    
                if (n2 <= specNSlopeBound(jNidx)) && (jNidx < 18)
                    jNidx = jNidx + 1;
                end
    
                if (n2 <= specNSlopeBound(jNidx)) && (jNidx >= 18)
                    jNidx = 18;
                end
    
                z0 = z2;     % n1 and z1 for next loop
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
    
        loudTDepChan(lStep) = N; % total loudness at current timeframe l
    end
    
    %%%% End of SQAT(/AARAE) code adapted
    specLoud(:, :, chan) = specLoudChan;
    loudTDep(:, chan) = loudTDepChan;
end  % end of for loop over channels

% power-averaged overall loudness
loudPowAvg = power(mean(loudTDep.^(1/log10(2)), 1), log10(2));

% 95th percentile overall loudness
loud5pcEx = prctile(loudTDep, 95, 1);

% power-averaged specific loudness as a function of Bark number
specLoudPowAvg = squeeze(power(mean(specLoud.^(1/log10(2)), 1), log10(2)));

% loudness level
switch adjustLoud
    case 'iso5321'
        model = 'zwicker';
    case 'iso5323'
        model = 'mgs';
    case 'ecma4182'
        model = 'sottek';
end
loudLevel = soneToPhon(loudTDep, model);
loudLvlPowAvg = soneToPhon(loudPowAvg, model);
loudLvl5pcEx = soneToPhon(loud5pcEx, model);

% time (s) corresponding with results output
timeOut = (0:(size(spectrL, 1) - 1))*timeStep;

%% Output assignment

loudness.loudTDep = loudTDep;
loudness.loudPowAvg = loudPowAvg;
loudness.loud5pcEx = loud5pcEx;
loudness.loudLevel = loudLevel;
loudness.loudLvlPowAvg = loudLvlPowAvg;
loudness.loudLvl5pcEx = loudLvl5pcEx;
loudness.specLoudPowAvg = specLoudPowAvg;
loudness.specLoud = specLoud;
loudness.barkAxis = barkAxis;
loudness.soundField = soundField;
loudness.adjustLoud = adjustLoud;
loudness.freqInMid = freqInMid;
loudness.freqInNom = freqInNom;
loudness.timeOut = timeOut.';

% used for development purposes
loudness.loudCorePowAvg = squeeze(power(mean(loudCore.^(1/log10(2)), 1), log10(2)));
loudness.loudCore = loudCore;

%% Output plotting

if outPlot
    % Plot figures
    % ------------
    dt = timeOut(2) - timeOut(1);
    freqAxis = bark2Hertz(barkAxis, 'volk');
    [~, ~, ~, fOctNom] = noctf(fLim, 1);
    [~, ~, ~, fOctNomAll] = noctf([25, 12600], 1);

    fDiffMat = abs(fOctNomAll(:) - fOctNom(:).');
    [~, minIdx] = min(fDiffMat, [], 1);
    
    fOctNomAllStr = ["31.5", "63", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"];
    fOctNomStr = fOctNomAllStr(minIdx);

    switch adjustLoud
        case 'iso5321'
            unitPlot = 'sone_{Z}';
            specUnitPlot =  strcat(unitPlot, '/Bark_{Z}');
        case 'iso5323'
            unitPlot = 'sone_{M-G-S}';
            specUnitPlot = strcat(unitPlot, '/Cam');
        case 'ecma4182'
            unitPlot = 'sone_{SHM}';
            specUnitPlot = strcat(unitPlot, '/Bark_{SHM}');
    end

    for chan = inChans:-1:1
        cmap_viridis = load('cmap_viridis.txt');
        % Plot results
        fig = figure;
        tiledlayout(fig, 2, 1);
        movegui(fig, 'center');
        ax1 = nexttile(1);
        surf(ax1, [timeOut, timeOut(end) + dt], freqAxis, [permute(specLoud(:, :, chan),...
                                              [2, 1, 3]), specLoud(end, :, chan).'],...
             'EdgeColor', 'none', 'FaceColor', 'interp');
        view(2);
        ax1.XLim = [timeOut(1), timeOut(end) + dt];
        ax1.YLim = [freqLower(1), freqUpper(end)];
        ax1.CLim = [0, ceil(max(specLoud(:, :, chan), [], 'all')*10)/10];
        ax1.YTick = fOctNom;
        ax1.YTickLabel = fOctNomStr;
        ax1.YScale = 'log';
        ax1.YLabel.String = 'Frequency, Hz';
        ax1.XLabel.String = 'Time, s';
        ax1.FontName =  'Arial';
        ax1.FontSize = 12;
        colormap(cmap_viridis);
        h = colorbar;
        set(get(h,'label'),'string', {'Specific Loudness,'; specUnitPlot});        
        
        ax2 = nexttile(2);
        plot(ax2, timeOut, loudPowAvg(1, chan)*ones(size(timeOut)), ':', 'color',...
             cmap_viridis(34, :), 'LineWidth', 1.5, 'DisplayName', "Power" + string(newline) + "time-avg");
        hold on
        plot(ax2, timeOut, loudTDep(:, chan), 'color', cmap_viridis(166, :),...
             'LineWidth', 0.75, 'DisplayName', "Time-" + string(newline) + "dependent");
        hold off
        ax2.XLim = [timeOut(1), timeOut(end) + (timeOut(2) - timeOut(1))];
        if max(loudTDep(:, chan)) > 0
            ax2.YLim = [0, 1.1*ceil(max(loudTDep(:, chan))*10)/10];
        end
        ax2.XLabel.String = 'Time, s';
        ax2.YLabel.String = strcat('Loudness,', " ", unitPlot);
        ax2.XGrid = 'on';
        ax2.YGrid = 'on';
        ax2.GridAlpha = 0.075;
        ax2.GridLineStyle = '--';
        ax2.GridLineWidth = 0.25;
        ax2.FontName = 'Arial';
        ax2.FontSize = 12;
        lgd = legend('Location', 'eastoutside', 'FontSize', 8);
        lgd.Title.String = "Overall";
    end  % end of for loop for plotting over channels
end  % end of if branch for plotting

end  % end of acousticQuasiLoudZwicker function