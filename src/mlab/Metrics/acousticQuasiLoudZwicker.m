function loudness = acousticQuasiLoudZwicker(spectrL, fLim, timeStep, axisF, soundField, adjustLoud, outPlot, recalibrate)
% loudness = acousticQuasiLoudZwicker(spectrL, fLim, timeStep, axisF,
%                                     soundField, adjustLoud, outPlot,
%                                     recalibrate)
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
% recalibrate : Boolean (default: false)
%   Indicates whether to apply a model-specific calibration to predict
%   absolute loudness values (derived from empirical data).
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
% Date last modified: 24/03/2026
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
        recalibrate (1, 1) {mustBeNumericOrLogical} = false
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
        calN = 1.213421616415088;
    case 'iso5323'
        calN = 1.293954956487481;
    case 'ecma4182'
        calN = 1.19668683901243;
end

% loudness transform threshold factor
threshFact = 0.25;

% loudness transform parameters (per critical band)
switch adjustLoud
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
end  % end of adjustLoud switch

% empirical calibration factors for complex
% sounds to recover absolute values closer to full
% psychoacoustic metrics
% NB: this recalibration means the reference sound will no longer equal 1 sone
if recalibrate
    switch adjustLoud
        case 'iso5321'
            calEmp = 0.82;  % this needs checking and updating

        case 'iso5323'
            calEmp = 1.0296789599477227;

        case 'ecma4182'
            calEmp = 0.6141829528096241;
    end
end

% Core loudness parameters for scaling adjustments
switch adjustLoud
    case 'iso5321'
        loudCoreScalLim = [0.175098518325143, 0.0998148899800973,...
                           0.0561937494710363, 0.0503003620878583,...
                           0.0323602219765719, 0.0214601301541275,...
                           0.0145579600580127, 0.0106902629505937,...
                           0.0120328553782146, 0.0156777674547530,...
                           0.0145268445060852, 0.00792251951842475,...
                           0.00530453058749494, 0.00715977949161123,...
                           0.0108711963376870, 0.0188389588811493,...
                           0.0438277330499417, 0.0692535984150476,...
                           0.0714553737495326, 0.0211134165088091;
                           27.8332159380776, 21.2463983253255,...
                           19.1459284746169, 18.3823397952096,...
                           17.5653691861714, 16.7503272771114,...
                           16.0584711550173, 15.5667156824197,...
                           15.7192882058850, 18.2757916953684,...
                           20.4705331364133, 18.3851539863524,...
                           17.4463354994553, 18.9373451514854,...
                           20.9017211077950, 23.3320877650296,...
                           27.2794008723165, 28.2581716202501,...
                           23.5149189707202, 10.2126762563267];

        polyCoeffs = [3.46154164787264e-06, 0.00115902620724796, 0.000908348090696937,...
                      0.00125059008443649, 0.000680464259160844, 0.000381612901601609,...
                      0.000250404904069992, 0.000307895392953073, 0.000370770653030356,...
                      0.000370770653030356, 0.000253806772389052, 0.000204153657667395,...
                      0.000116912083670013, 0.000124352541267608, 0.000166448922257427,...
                      0.000484847673433429, 0.000879751536044967, 0.000909345764564276,...
                      0.000522148404621712, 0.000208797319808608;
                      -0.00114792490037041, -0.00484549495012134, -0.00346106937185142,...
                      -0.00306582238333499, -0.00255740228537038, -0.00245083618521589,...
                      -0.00215740686274558, -0.00184701927288191, -0.00189366867382637,...
                      -0.00189366867382637, -0.00184025143910359, -0.00157341566150419,...
                      -0.00164106124813877, -0.00189288456298325, -0.00220988694211908,...
                      -0.00231272681598800, -0.00252227743347211, -0.00163103198459638,...
                      0.000725192584004455, 0.00120141808594060;
                      -0.0107572259682556, -0.0216831709605858, -0.0199503065136557,...
                      -0.0322758543645889, -0.0202998722763244, -0.0126333392687931,...
                      -0.00933153105471622, -0.0115007248697960, -0.0126058004391672,...
                      -0.0126058004391672, -0.00913000773986548, -0.00940237723359883,...
                      -0.00786739562303758, -0.00773095294646210, -0.00863282700635499,...
                      -0.0170487449794070, -0.0261267286816691, -0.0281851220710378,...
                      -0.0238906686292304, -0.00723592624190622;
                      0.0960537985518171, 0.118672875233916, 0.111336651243432,...
                      0.130133891464595, 0.114874037278364, 0.105238379315609,...
                      0.0933133601822083, 0.0893081868136628, 0.0863416929336673,...
                      0.0863416929336673, 0.0801290638010130, 0.0789428286903867,...
                      0.0819485300104738, 0.0837324848247202, 0.0913019376436802,...
                      0.0993924998215976, 0.0987350870968390, 0.0781892225672751,...
                      0.0448080957714798, 0.00343201959298394;
                      1.13398696438849, 1.12472183680387, 1.11006617352212,...
                      1.12267702135544, 1.10766334397886, 1.09961670412122,...
                      1.10837503442643, 1.12474743102009, 1.13255905488564,...
                      1.13255905488564, 1.13070313922545, 1.13161344455377,...
                      1.12325377295784, 1.12256168883824, 1.11750024738604,...
                      1.13179880735164, 1.17314654781390, 1.20975095639261,...
                      1.22642463134023, 1.22679079317417];

        xShift = 0;  % null shift for ISO 532-1 log polynomial fit

    case 'iso5323'
        loudCoreScalLim = [0, 0.00530830737662207, 0.0192031892338177,...
                           0.0209111591599521, 0.0128984964763127,...
                           0.0104915273945527, 0.0116827917322032,...
                           0.00811613906542941, 0.00730727308003143,...
                           0.00636456516716606, 0.00954914811672541,...
                           0.0102371893640083, 0.0100504773914230,...
                           0.0105945560510230, 0.0125674814276614,...
                           0.0125443298326316, 0.0106942015417099,...
                           0.00841349456970890, 0.0134021791395764, 0.0176055747172162;
                           14.5713717505921, 20.0186207046044,...
                           24.3422225721081, 23.1537750300632,...
                           25.7105391976310, 25.5632954402285,...
                           27.7329515398835, 27.9821789123069,...
                           27.6712025848049, 30.0587852761365,...
                           37.7528978299404, 36.9301593441850,...
                           36.9195806911093, 37.4360689086251,...
                           38.5028501361481, 35.8850259099337,...
                           29.0401693817834, 22.3127471294020,...
                           17.5043760969581, 9.60794408548125];

        polyCoeffs = [-0.0315555801173618, -0.0235155111682275, -0.00876096637958047,...
                      -0.00461946156700320, -0.0128163078458873, -0.0174003747974420,...
                      -0.0144126636878137, -0.0144819651160435, -0.0176099844739473,...
                      -0.0176099844739473, -0.0101402738663035, -0.0100053187922727,...
                      -0.0117975068960418, -0.0114079876138620, -0.00987023464425238,...
                      -0.00801237581713489, -0.0167425135148038, -0.0286038654484740,...
                      -0.0275718886697226, -0.0446111736908394;
                      0.272898875453994, 0.226141813402174, 0.100783539800245,...
                      0.0585363180759330, 0.146117398400177, 0.187896147115585,...
                      0.159720565928761, 0.162775919044573, 0.201559150313230,...
                      0.201559150313230, 0.128959100607096, 0.123341273427516,...
                      0.141468320936720, 0.137233388104629, 0.123225173355870,...
                      0.103599501365731, 0.196253756184640, 0.310778890388836,...
                      0.271813841877534, 0.332446752065889;
                      -0.778006230398953, -0.740575990117963, -0.393336181560437,...
                      -0.262781629784085, -0.565974640296568, -0.687355110768661,...
                      -0.599826125785208, -0.627029241220489, -0.782991128640219,...
                      -0.782991128640219, -0.555142814062536, -0.518530814708928,...
                      -0.581639697674929, -0.572116003342879, -0.540246178581444,...
                      -0.488705882137998, -0.817725947432973, -1.17754879591352,...
                      -0.929504774083775, -0.810014921346866;
                      1.09310213153993, 1.21690500206716, 0.914064764243455,...
                      0.793625036782538, 1.17684011028206, 1.30342601070191,...
                      1.19400644114233, 1.26692294231247, 1.48590109031404,...
                      1.48590109031404, 1.23247301731577, 1.15619750543038,...
                      1.24244915221045, 1.23671441009310, 1.22026888327745,...
                      1.20927878638716, 1.59483382755989, 1.94535723221769,...
                      1.41122356398346, 0.821846998291404;
                      0.142170260463131, 0.374853893016304, 0.645503513948269,...
                      0.410437564534149, 0.389933440723440, 0.419750124030468,...
                      0.544545848902599, 0.482739931140136, 0.474773701207036,...
                      0.474773701207036, 0.609416397321254, 0.588748217131579,...
                      0.506907680129803, 0.561059682577414, 0.554852082306551,...
                      0.472847255737247, 0.549210638161651, 0.425449811688886,...
                      0.367695270514220, 0.480884272265022];

        xShift = 1;  % shift for ISO 532-3 log polynomial fit

    case 'ecma4182'
        loudCoreScalLim = [0, 0.0146423830980930, 0.0401421283871121,...
                           0.0264452523783859, 0.0260063706293114,...
                           0.0270032954442271, 0.0287488839854639,...
                           0.0305306502376336, 0.0306362496450788,...
                           0.0368871112603854, 0.0318958841234506,...
                           0.0220432101869425, 0.0139813205583574,...
                           0.0159376609929878, 0.0177891820900005,...
                           0.0231884383366289, 0.0266651346569812,...
                           0.0311287979457835, 0.0394127546716581, 0.0291541517395267;
                           5.12035926747327, 6.72365689831412,...
                           7.92283899060251, 7.74024319015972,...
                           8.32590649318051, 8.70378353188929,...
                           9.25191870233799, 9.25649634953904,...
                           8.84900560365458, 8.57309331460825,...
                           9.81329560391793, 9.31097256099634,...
                           8.97029603765363, 9.49786753160517,...
                           9.62114709441200, 10.2322084436148,...
                           10.0256382511768, 8.82618293558970,...
                           7.73558515861310, 4.40280304260254];

        polyCoeffs = [-0.0736891438323439, -0.0159530506355199, 0.0212550722617645,...
                      -0.0200209140249662, -0.0226540435091588, -0.0395070756417713,...
                      -0.0303603139392496, -0.0240070928379679, 0.00131842305195210,...
                      0.00131842305195210, 0.0184672510479861, -0.0476845602661092,...
                      -0.0569338047946894, -0.0308625608385108, -0.0718774569482032,...
                      -0.0555031470988601, -0.0269592580328595, -0.0409248088421269,...
                      -0.0407963332912421, -0.0389642975658443;
                      0.530698179354883, 0.141627611044728, -0.152474400844396,...
                      0.168121212537712, 0.199606979004556, 0.320223322274041,...
                      0.260042548319418, 0.209242752348185, 0.0197713523776238,...
                      0.0197713523776238, -0.113303575398161, 0.382379630581085,...
                      0.470367290019549, 0.270066922632544, 0.515621465942002,...
                      0.404038059625160, 0.234919582328986, 0.323611017174446,...
                      0.283655227677970, 0.113449591729784;
                      -1.30022426384128, -0.430137862993729, 0.329666764298806,...
                      -0.537085334793378, -0.633459935768959, -0.912381421553071,...
                      -0.798347250285430, -0.668931362568729, -0.220595747675774,...
                      -0.220595747675774, 0.123716962668953, -1.09617281773497,...
                      -1.38230500354192, -0.888962952125888, -1.32460375705285,...
                      -1.07990833150335, -0.806930522071779, -1.00741195504106,...
                      -0.783006683287310, -0.0275699333975544;
                      1.44055337832766, 0.784650198067210, 0.125789648094553,...
                      1.05518882201568, 1.15263650146115, 1.38370167270274,...
                      1.33902047956841, 1.22274940362274, 0.883083282911430,...
                      0.883083282911430, 0.526466047408033, 1.64210224054734,...
                      1.99045255900662, 1.57051418723404, 1.80564205150661,...
                      1.62465493650303, 1.52666192816348, 1.66812217236099,...
                      1.23011075454101, 0.140621133557864;
                      0.0953333449232013, 0.444241924569629, 0.801522281515868,...
                      0.331072353747749, 0.395343001991298, 0.414240865116689,...
                      0.486045758592103, 0.529280508893734, 0.776069205096691,...
                      0.776069205096691, 0.950444998482004, 0.533943294885853,...
                      0.391928101633608, 0.494126668846407, 0.414777893543282,...
                      0.452220605922443, 0.616625097506653, 0.573531007383365,...
                      0.476739060897969, 0.629828800473469];

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
        z = 0.1;
        if CB1
            z1 = 0;  % critical band rate starts at 0
        else
            z1 = min(barkNRange);
        end
        iz = 1;  % critical band index
        jNidx = 18;  % index of loudness for slope calculation (Table A.9)
    
        for iBand = 1:barkCount + 1  % specific loudness
    
            % Determines where to start on the slope
            ig = iBand + barkI1 - 2;  % the addition accounts for a start band above band 1
    
            % steepness of upper slope (specNSlope) for bands above 8th one are identical
            if ig > 8
                ig = 8;
            end
    
            while z1 < barkNRange(iBand)
    
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
                    N = N + n2*(z2 - z1);
                    zk = z;                     % initialisation of zk
    
                    while (zk <= z2)
                        specLoudChan(lStep, iz) = n2;
                        iz = iz + 1;
                        zk = zk + 0.1;
                    end
    
                    z = zk;
    
                else  % if slope loudness > core loudness
                    
                    % decision whether the critical band in question is completely
                    % or partly masked by accessory loudness
    
                    n2 = specNSlopeBound(jNidx);
    
                    if n2 < loudCore1(lStep, iBand)
                        n2 = loudCore1(lStep, iBand);
                    end
    
                    dz = (n1 - n2) / specNSlope(jNidx, ig);
                    z2 = z1 + dz;
    
                    if z2 > barkNRange(iBand)
                        z2 = barkNRange(iBand);
                        dz = z2 - z1;
                        n2 = n1 - dz*specNSlope(jNidx, ig);
                    end
    
                    N = N + dz*(n1 + n2)/2;
                    zk = z;                     % initialisation of zk
    
                    while (zk <= z2)
                        specLoudChan(lStep, iz) = n1 - (zk - z1)*specNSlope(jNidx, ig);
                        iz = iz + 1;
                        zk = zk + (1/10);
                    end
    
                    z = zk;
    
                end
    
                if (n2 <= specNSlopeBound(jNidx)) && (jNidx < 18)
                    jNidx = jNidx + 1;
                end
    
                if (n2 <= specNSlopeBound(jNidx)) && (jNidx >= 18)
                    jNidx = 18;
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
    
        loudTDepChan(lStep) = N; % total loudness at current timeframe l
    end
    
    %%%% End of SQAT(/AARAE) code adapted
    specLoud(:, :, chan) = specLoudChan;
    loudTDep(:, chan) = loudTDepChan;
end  % end of for loop over channels

% apply calibrations, if using
if recalibrate
    specLoud = calEmp*specLoud;
    loudTDep = calEmp*loudTDep;
end

% power-averaged overall loudness
loudPowAvg = power(mean(loudTDep.^(1/log10(2)), 1), log10(2));

% 95th percentile overall loudness
loud5pcEx = prctile(loudTDep, 95, 1);

% power-averaged specific loudness as a function of Bark number
specLoudPowAvg = squeeze(power(mean(specLoud.^(1/log10(2)), 1), log10(2)));

% loudness level
if recalibrate && ~strcmp(adjustLoud, 'iso5321')
    switch adjustLoud
        case 'iso5323'
            model = 'mgs';
        case 'ecma4182'
            model = 'sottek';
    end
    loudLevel = soneToPhon(loudTDep, model);
    loudLvlPowAvg = soneToPhon(loudPowAvg, model);
    loudLvl5pcEx = soneToPhon(loud5pcEx, model);
else
    loudLevel = soneToPhon(loudTDep, 'zwicker');
    loudLvlPowAvg = soneToPhon(loudPowAvg, 'zwicker');
    loudLvl5pcEx = soneToPhon(loud5pcEx, 'zwicker');
end

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