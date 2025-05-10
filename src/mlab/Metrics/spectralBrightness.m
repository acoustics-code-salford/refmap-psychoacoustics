function brightness = spectralBrightness(spectr, octN, fLim, axisT, tStep, method, outPlot)
% brightness = spectralBrightness(spectr, octN, fLim, axisT, tStep, method, outPlot)
%
% Returns spectral brightness (approximation to psychoacoustic sharpness)
% using fractional octave band spectral input, as a single unweighted
% spectrum, or a matrix of (e.g., time-dependent) spectra.
%
% Inputs
% ------
% spectr : vector or 2D matrix
%          the continguous input spectrum or spectra for processing
%
% octN : integer (default: 3)
%        the octave fraction denominator (e.g., 3 indicates 1/3-octave
%        spectra)
%
% fLim : vector (default: [20, 20e3])
%        the frequency limits for the input spectra (these are
%        automatically matched to nearest exact band centre-frequencies for
%        the selected octN resolution
%
% axisT : integer (default: 1)
%         the axis that is not the frequency band axis (e.g., the time
%         axis)
%
% tStep : float (default: 0.5 s)
%         the time step (seconds) for time-dependent spectra. This is only
%         used for plotting purposes.
%
% method : keyword string (default: 'vonbismarck')
%          the sharpness weighting function to approximate (options
%          comprise 'vonbismarck', 'aures' or 'widmann'.
%
% outPlot : Boolean (default: false)
%           determines whether to plot outputs from the calculations
% 
% Returns
% -------
% brightness : structure
%              contains the output
%
% brightness contains the following outputs:
%
% dBSpecTarget : matrix
%                sound pressure level spectrogram for the input target
%                signal, with dimensions [timeOut, freqBands, targetChans]
%
% dBSpecMasker : matrix
%                sound pressure level spectrogram for the input masker
%                signal, with dimensions [timeOut, freqBands, maskerChans] 
%
% dBSpecDiscTarget : matrix
%                    sound pressure level spectrogram for the input target
%                    signal detectability-discounted, with dimensions
%                    [timeOut, freqBands, targetChans]
%
% dBATDepTarget : matrix or vector
%                 time-dependent A-weighted sound pressure level for the
%                 input target signal, with dimensions [timeOut, targetChans]
%
% dBATDepTarget : matrix or vector
%                 time-dependent A-weighted sound pressure level for the
%                 input masker signal, with dimensions [timeOut, targetChans]
%
% dBATDepDiscTarget : matrix or vector
%                     time-dependent A-weighted sound pressure level for the
%                     input target signal detectability-discounted, with
%                     dimensions [timeOut, targetChans]
%
% dBATDepDiscount : vector
%                   time-dependent detectability discount values for the
%                   A-weighted sound pressure levels of target signal vs
%                   masker signal, with dimensions [timeOut, targetChans]
%
% LAETarget : vector
%             overall A-weighted sound exposure level for each input target
%             signal channel
%
% LAEMasker : vector
%             overall A-weighted sound exposure level for each input masker
%             signal channel
%
% LAEDiscTarget : vector
%                 overall detectability-discounted A-weighted sound
%                 exposure level for each input target signal channel
%
% LAeqTarget : vector
%              overall A-weighted time-averaged sound level for each input
%              target signal channel
%
% LAeqMasker : vector
%              overall A-weighted time-averaged sound level for each input
%              masker signal channel
%
% LAeqDiscTarget : vector
%                  overall detectability-discounted A-weighted
%                  time-averaged sound level for each input target signal
%                  channel
%
% dBADiscount : vector
%               overall detectability discount for A-weighted levels for
%               each input target signal channel
%
% detectabilitydB : matrix
%                   detectability spectrogram for the input target signal
%                   vs masker signal (i.e., 10log_10[d']), with dimensions
%                   [timeOut, freqBands, targetChans]
%
% detectTDepMaxdB : matrix or vector
%                   band-maximum time-dependent detectability dB, with
%                   dimensions [timeOut, targetChans]
%
% detectTDepIntdB : matrix or vector
%                   band-integrated time-dependent detectability dB, with
%                   dimensions [timeOut, targetChans]
%
% detectMaxdB : vector
%               overall maximum detectability, dB, for each input target
%               signal channel
%
% detectIntdB : vector
%               overall intergated detectability, dB, for each input target
%               signal channel
%
% detectMaxPcdB : structure
%                 contains percent-based aggregated metrics for maximum
%                 detectability, comprising:
% Ex50 : vector
%        50% exceeded (median) dB, for each input target signal channel
%
% detectIntPcdB : structure
%                 contains percent-based aggregated metrics for maximum
%                 detectability, comprising:
% Ex50 : vector
%        50% exceeded (median) dB, for each input target signal channel
%
% freqBands : vector
%             1/3-octave band centre-frequencies for input freqBandRange
%
% timeOut : vector
%           the window centre times for the spectrograms
%
% Assumptions
% -----------
% The input is a fractional octave band unweighted energy time-averaged
% sound level (LZeq) spectrum or series of (time-averaged) LZeq spectra in
% dB re 2e-5 Pa.
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
% Date created: 22/04/2025
% Date last modified: 22/04/2025
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
        spectr (:, :) double {mustBeReal}
        octN (1, 1) double {mustBePositive, mustBeInteger} = 3
        fLim (1, 2) double {mustBeInRange(fLim, 19, 20000)} = [20, 20000]
        axisT (1, 1) {mustBeInteger, mustBeInRange(axisT, 1, 2)} = 1
        tStep (1, 1) double {mustBePositive} = 0.5
        method (1, :) string {mustBeMember(method,...
                                           {'aures',...
                                            'vonbismarck',...
                                            'widmann'})} = 'vonbismarck'
        outPlot {mustBeNumericOrLogical} = false
    end

%% Load path
addpath(genpath(fullfile("refmap-psychoacoustics", "src", "mlab")))

%% Input checks
if axisT == 2
    spectr = spectr.';
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

% estimated extension of alpha for 16 and 20 kHz bands
alphaExt = [alpha, 0.387, 0.480];

% ISO 226:2023 Lu values
Lu = [-31.5, -27.2, -23.1, -19.3, -16.1, -13.1, -10.4, -8.20, -6.30,...
      -4.60, -3.20, -2.10, -1.20, -0.50, 0.00, 0.40, 0.50, 0.00, -2.70,...
      -4.20, -1.20, 1.40, 2.30, 1.00, -2.30, -7.20, -11.2, -10.9, -3.50];

% estimated extension of Lu for 16 and 20 kHz bands
LuExt = [Lu, -2, -7];

% normal hearing from ISO 226:2023 (20 Hz - 12.5 kHz)
hearThresholds226 = [78.1 68.7 59.5 51.1 44.0 37.5 31.5 26.5 22.1...
                     17.9 14.4 11.4 8.6 6.2 4.4 3.0 2.2 2.4 3.5...
                     1.7 -1.3 -4.2 -6.0 -5.4 -1.5 6.0 12.6 13.9 12.3];

% free-field frontal tone hearing thresholds for 18-25 year-olds with
% normal hearing from ISO 389-7:2019 (16 - 18 kHz)
hearThresholds3897 = [40.2, 70.4];

% free-field frontal tone hearing thresholds for 18-25 year-olds with
% normal hearing from ISO 226:2023 | ISO 389-7:2019 (20 - 18 kHz), with 20
% kHz band estimated from figure 1 in ISO 389-7:2019
hearThresholdsFFTone = [hearThresholds226, hearThresholds3897];

% diffuse field narrowband noise hearing thresholds for 18-25 year-olds with
% normal hearing from ISO 389-7:2019 (20 Hz - 16 kHz)
% NOTE: watch out for non-standard frequencies in the standard Table 1!
hearThresholdsDF3897 = [78.1, 68.7, 59.5, 51.1, 44.0, 37.5, 31.5, 26.5, 22.1,...
                        17.9, 14.4, 11.4, 8.4, 5.8, 3.8, 2.1, 1.0, 0.8, 1.9,...
                        0.5, -1.5, -3.1, -4.0, -3.8, -1.8, 2.5, 6.8, 8.4, 14.4,...
                        43.7];

% estimated full range diffuse field narrowband noise hearing thresholds
% for 18-25 year-olds with normal hearing 
hearThresholdsDFNoise = [hearThresholdsDF3897, 70.4];

%% Weighting function
% fractional octave frequencies
[fm, f1, f2, fn] = noctf(fLim, octN);

% fractional octave bandwidths
bandwidthfOct = f2 - f1;

% interpolate functions over all frequencies
[fm3, ~, ~, ~] = noctf([20 20e3], 3);
if octN == 3
    alphaUse = alphaExt;
    LuUse = LuExt;
    hTUse = hearThresholdsDFNoise;
else
    alphaFit = fit(fm3.', alphaExt.', 'cubicspline');
    LuFit = fit(fm3.', LuExt.', 'cubicspline');
    hTFit = fit(fm3.', hearThresholdsDFNoise, 'cubicspline');
    alphaUse = alphaFit(fm);
    LuUse = LuFit(fm);
    hTUse = hTFit(fm);
end

alpha1k = repmat(alphaUse(fm==1e3), 1, size(alphaUse, 2));
Lu1k = repmat(LuUse(fm==1e3), 1, size(LuUse, 2));
hT1k = repmat(hTUse(fm==1e3), 1, size(hTUse, 2));

%% Signal processing
soundLevels = loudLfISO226_2023(alphaUse, LuUse, hTUse, (20:10:90).');
weighting = -(soundLevels - loudLfISO226_2023(alphaUse(fm==1e3), LuUse(fm==1e3), hTUse(fm==1e3), (20:10:90).'));

end  % end of function

