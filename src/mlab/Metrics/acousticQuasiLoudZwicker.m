function loudness = acousticQuasiLoudZwicker(spectrL, fLim, axisF, soundField, adjustLoud)
% loudness = acousticQuasiLoudZwicker(spectrL, fLim, axisF, soundField,
%                                      adjustLoud)
%
% Returns quasi-loudness using spectral elements of ISO 532-1 Zwicker
% loudness model for arbitrary spectra. The temporal weighting of the
% model is disregarded. The input spectra must be 1/3 octave band
% resolution Leq, limited to the 25 Hz-12.5 kHz range.
%
% Optional modifications available comprise an adjustment to the spectral
% levels to improve agreement with the 2023 ISO 226 equal loudness
% contours ('iso226'), or (alternatively) the application of the
% ECMA-418-2:2025 outer-middle ear filter responses (in 1/3-octaves)
% instead of the ISO 532-1:2017 outer ear transmission, and an
% approximated non-linear transformation following ECMA-418-2:2025
% ('ecma4182'), instead of Zwicker's original.
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
% adjustLoud : keyword string (default: 'none')
%   Indicates whether to apply adjustments for:
%   ('iso226') the differences between 1987 ISO 226 equal-loudness contours
%   (which ISO 532-1 models) and 2023 ISO 226 equal-loudness contours, or;
%   ('ecma4182') the outer-middle ear filter response from ECMA-418-2:2024
%   omitting the ISO 532-1:2017 critical band ear transmission a0, and
%   adapting the loudness transformation to agree more closely with the
%   ECMA-418-2:2024 transformation.
%   The default option ('none') applies no adjustment, so follows ISO 532-1
%   more closely (for closer agreement with Zwicker's model).
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
% Date last modified: 27/11/2025
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
                                               {'none',...
                                                'iso226',...
                                                'ecma4182'})} = 'none'
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

% ISO 532-1 Table A.6: critical band level at threshold in quiet
levelThresQ = [30, 18, 12, 8, 7, 6, 5, 4,...
               3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3];

% ISO 532-1 Table A.7: critical bandwidth vs. 1/3-octave bandwidth
% adjustment
critBWAdjust = [-0.25, -0.6, -0.8, -0.8, -0.5, 0, 0.5, 1.1, 1.5, 1.7,...
                1.8, 1.8, 1.7, 1.6, 1.4, 1.2, 0.8, 0.5, 0, -0.5];

if strcmp(adjustLoud, 'ecma4182')
    % ISO 523-2 Table 2
    levelThresQTOB = [28.18, 23.9, 19.2, 15.68, 12.67, 10.09, 8.08, 6.3,...
                      5.3, 4.5, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63,...
                      3.63, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63, 3.63];

    levelThresQCB = [min(levelThresQTOB(1:3)), min(levelThresQTOB(4:6)),...
                     min(levelThresQTOB(7:8)), levelThresQTOB(9:end)];

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

% level difference between 2023 ISO 226 equal loudness contours (10-100
% phon, 25-12500 Hz) and 1987 ISO 226 contours
d2023v1987 = [ 8.4, 10.4, 11.0, 11.0, 10.6,  9.8,  8.5, 6.9, 4.8,  0.0;
                 9.5, 12.0, 12.8, 12.8, 12.4, 11.5, 10.1, 8.4, 6.3,  3.6;
                10.2, 12.8, 13.7, 13.7, 13.1, 12.1, 10.6, 8.8, 6.6,  4.0;
                10.4, 13.1, 13.9, 13.8, 13.0, 11.9, 10.4, 8.5, 6.2,  3.5;
                10.3, 13.1, 13.8, 13.5, 12.7, 11.5,  9.9, 7.9, 5.6,  2.9;
                10.0, 12.9, 13.5, 13.3, 12.5, 11.3,  9.7, 7.7, 5.4,  2.7;
                 9.5, 12.4, 13.2, 13.1, 12.3, 11.2,  9.7, 7.8, 5.5,  2.7;
                 9.1, 12.0, 12.9, 12.9, 12.3, 11.2,  9.8, 7.9, 5.6,  2.9;
                 8.2, 11.1, 12.2, 12.3, 12.0, 11.0,  9.7, 7.9, 5.7,  3.1;
                 7.0,  9.9, 11.2, 11.6, 11.3, 10.5,  9.4, 7.8, 5.9,  3.3;
                 5.9,  8.8, 10.0, 10.6, 10.5,  9.9,  9.0, 7.5, 5.7,  3.3;
                 4.9,  7.5,  8.8,  9.4,  9.4,  9.1,  8.2, 6.9, 5.3,  3.2;
                 3.5,  5.8,  7.1,  7.7,  7.8,  7.5,  6.8, 5.9, 4.5,  2.6;
                 2.3,  4.2,  5.3,  5.8,  6.0,  5.9,  5.5, 4.7, 3.6,  2.2;
                 1.2,  2.5,  3.3,  3.7,  3.9,  3.8,  3.4, 3.0, 2.2,  1.3;
                 0.1,  0.7,  1.0,  1.3,  1.4,  1.4,  1.3, 1.1, 0.8,  0.4;
                 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0;
                 1.5,  1.5,  1.5,  1.6,  1.8,  2.0,  2.3, 2.7, 3.1,  3.7;
                 1.8,  2.3,  2.6,  3.0,  3.3,  3.7,  4.2, 4.9, 5.7,  6.7;
                 0.1,  0.5,  0.7,  1.1,  1.6,  2.3,  3.1, 4.0, 5.2,  6.5;
                -0.5,  0.1,  0.4,  0.9,  1.5,  2.3,  3.3, 4.4, 5.8,  7.3;
                 0.4,  1.2,  1.7,  2.4,  3.0,  3.9,  5.0, 6.2, 7.5,  9.1;
                 1.7,  2.8,  3.5,  4.2,  4.9,  5.9,  6.8, 8.0, 9.2, 10.7;
                 2.8,  3.7,  4.5,  5.2,  5.8,  6.6,  7.4, 8.3, 9.2, 10.3;
                 2.2,  3.0,  3.8,  4.3,  4.9,  5.5,  5.9, 6.5, 7.0,  7.4;
                 0.5,  2.0,  3.1,  3.8,  4.4,  4.8,  5.0, 5.0, 4.9,  4.5;
                 2.7,  5.5,  7.0,  7.9,  8.3,  8.2,  7.8, 7.0, 5.7,  0.0;
                 6.9, 10.3, 12.0, 12.7, 12.5, 11.6,  9.7, 7.0, 3.2,  0.0];

phonRange = 10:10:100;

% Outer and middle ear filter functions corresponding with ECMA-418-2:2024,
% converted to 1/3-octave bands (25 - 12500 Hz)
outMidFree = [-22.31, -20.32, -18.45, -16.52, -14.64, -12.77, -11.17,...
               -9.62, -8.16, -6.95, -5.74, -4.39, -2.89, -1.31, 0.19,...
                1.27, 0.00, -2.71, -1.72, 1.31, 3.55, 5.24, 6.17,...
                3.78, -2.5, -8.34, -10.08, -10.64];

outMidDiffuse = [-22.31, -20.32, -18.45, -16.52, -14.64, -12.77, -11.18,...
                 -9.62, -8.17, -6.98, -5.79, -4.48, -3.09, -1.74, -0.76,...
                 -0.17, 0.24, 1.03, 2.29, 3.57, 4.72, 5.87, 6.52, 4.00,...
                 -2.37, -8.27, -10.05, -10.62];

midOnly = [-22.31, -20.32, -18.45, -16.52, -14.64, -12.77, -11.18,...
            -9.62, -8.17, -6.97, -5.77, -4.44, -3.00, -1.55, -0.4,...
             0.37, 0.63, 0.42, -0.17, -1.03, -2.03, -3.05, -4.00,...
            -4.77, -5.32, -5.66, -6.01, -7.42];

% Calibration constants (to ensure 1 sone corresponds with 1 kHz sinusoid at
% 40 dB in free-field)
switch adjustLoud
    case 'ecma4182'
        calN = 0.915728923552739;
    case 'iso226'
        calN = 1.209768056417869;
    case 'none'
        calN = 1.192929376089246;
end

% loudness transform threshold factor
threshFact = 0.25;

% loudness transform parameters per critical band for ecma4182 option
if strcmp(adjustLoud, 'ecma4182')

    tFormParams = [0.0216211014619521, 0.244893066697970, 0.249652468585476,...
                   0.132769990262922, 0.144338387038054, 0.138547904063839,...
                   0.127851720890203, 0.122988446348321, 0.109432262457469,...
                   0.140110970198732, 0.127532531503563, 0.128824953232883,...
                   0.125423326764015, 0.113713974586521, 0.114221022126836,...
                   0.115674103067031, 0.117045217317147, 0.123614083758182,...
                   0.128292534021925, 0.129388081356235;
                   12.9997784114591, 0.189689480014467, 27.5041570037957,...
                   14.0747011015335, 26.9618721060537, 14.4981572709487,...
                   3.96076147304715, 25.0327217497616, 1.16375666073087,...
                   21.7640099113582, 15.7816186599824, 1.62680264429789,...
                   2.11564296029020, 19.7133969292022, 12.2975292033567,...
                   14.9609171919786, -2.32510656585183, 10.4449368579285,...
                   11.5131700928274, 10.9776844794432;
                   0.550004650415479, 0.100080757265366, 0.346122712797352,...
                   0.550189114019556, 0.477856625731104, 0.456954428002810,...
                   0.374757323668428, 0.605408817611497, 0.216031615292129,...
                   0.582191786753812, 0.557145246877237, 0.342554465472874,...
                   0.370810435822780, 0.517159004220563, 0.831437862495454,...
                   0.826655742289538, 0.846147528058533, 0.813828381081569,...
                   0.769281361498644, 0.798452453256665;
                   0.234168807782939, 0.142865229874899, 0.147443484602303,...
                   0.191253199663201, 0.180718235313210, 0.182356695592064,...
                   0.187980849430008, 0.189686205126263, 0.210018882817866,...
                   0.202658092740540, 0.208829051454368, 0.268132695041064,...
                   0.203618421505789, 0.193542579976283, 0.187146691691566,...
                   0.185856331078352, 0.187151425242945, 0.187817210563186,...
                   0.184885920516634, 0.174287307872786;
                  -0.0314967447291350, -0.0506714665075209, -0.0376645433019592,...
                  -0.0351674378501139, -0.0194388857021200, -0.0152717939057543,...
                  -0.0134450413035557, -0.0116995448707150, -0.0213778718674300,...
                  -0.0348463270274590, -0.0347986379310834, -0.0925973329499940,...
                  -0.0247469217021565, -0.00833770842655125, -0.00310967338030957,...
                  -0.00373157927397809, -0.00718325628212586, -0.0137749197165944,...
                  -0.0165767904749479, -0.0119545792193787;
                   35.0336084516183, 83.3368871957613, 84.1501555884328,...
                   100.218714680987, 88.5903082350575, 87.5515649841792,...
                   90.4705999145900, 91.9220231584999, 118.982821284743,...
                   124.490670419392, 126.031894509287, 155.847920405969,...
                   119.855891847784, 100.362239143157, 93.0345292038482,...
                   93.0241169457276, 100.002846880705, 104.574081521952,...
                   103.922064899048	97.2030112052792;
                  -25.5529483888763, -47.2052688039484, -40.6420481056860,...
                  -57.8551311301582, -31.8981349724947, -29.2458884789345,...
                  -26.7512699340362, -26.7056271886970, -33.8428265154256,...
                  -49.1980895148980, -45.8712359428293, -42.3574478012233,...
                  -26.9826496760918, -16.4844827430767, -6.36146379509605,...
                  -1.82053441691555, -9.17083504578188, -15.2503110339556,...
                  -16.5317770116393, -11.7668051136539];

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
    rangePhons = length(phonRange);
    
    % loop though 1/3-oct bands
    for kk = fmAllI1:fmAllI2
    
        % only proceed if any band weightings are non-zero
        if any(lowFWeights(kk, :)) && ~strcmp(adjustLoud, 'ecma4182')
    
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

        % if adjusting for difference between 1987 vs 2023 equal loudness
        % contours
        if strcmp(adjustLoud, 'iso226')
            levelWeightedRep = repmat(levelWeighted(:, kk), 1, rangePhons);
    
            % check which levels meet corresponding criterion
            meetCriterion = levelWeightedRep <= phonRange;
    
            % get column indices for lowest criterion match
            [~, dCriterionMatch] = max(meetCriterion == 1, [], 2);
    
            % repeat differences matrix for each time step
            d2023v1987Rep = repmat(d2023v1987(kk, :), nTimeSteps, 1);
    
            % adjust levels according to matched criteria
            levelWeighted(:, kk) = levelWeighted(:, kk)...
                                   - d2023v1987Rep(sub2ind(size(d2023v1987Rep),...
                                   (1:nTimeSteps).', dCriterionMatch));
    
        end  % end of if branch for 1987 vs 2023 equal loudness adjustment
    
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
    excitCBands = excitation - critBWAdjust(barkI1:barkI2);
    maskLTQ = excitation > levelThresQ(barkI1:barkI2);
    levelThresQRep = repmat(levelThresQ(barkI1:barkI2), size(excitation, 1), 1);

    if strcmp(adjustLoud, 'ecma4182')
        
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

    else
        loudCore = calN.*((1 - threshFact...
                           + threshFact.*10.^((excitCBands...
                                               - levelThresQRep)/10)).^0.25...
                          - 1).*0.0635.*(10.^(0.1*levelThresQRep)).^0.25;
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