function detectDiscount = acousticDetectDiscount(signalTarget, sampleRateTarget, signalMasker,...
                                                 sampleRateMasker, axisTarget, axisMasker,...
                                                 timeStep, freqBandRange, outPlot)
% detectDiscount = acousticDetection(signalTarget, sampleRateTarget, axisTarget,
%                                    signalMasker, sampleRateMasker, axisMasker)
%
% Returns detectability and discounted sound levels from input target
% source and masker signals based on the detectability model originally
% developed by Bolt, Beranek and Newman consulting engineers, and developed
% further by NASA, with discounting of target source sound levels using the
% technique developed by NASA.
%
% Fidell, S et al, 1974. Prediction of aural detectability of noise
% signals, Human Factors 16(4), 373-383.
% https://doi.org/10.1177/001872087401600405
%
% Sneddon, M et al, 2003. Laboratory study of the noticeability and
% annoyance of low signal-to-noise ratio sounds, Noise Control Engineering
% Journal, 51(5), 300-305.
% https://doi.org/10.3397/1.2839726
%
% Christian, A, 2021. A construction for the prediction of noise-induced
% annoyance in the presence of auditory masking, 181st Meeting of the
% Acoustical Society of America.
% https://ntrs.nasa.gov/citations/20210024824
%
% Rizzi, SA et al, 2024.
%
%
% Inputs
% ------
% signalTarget : vector or 2D matrix
%                the input target signal as single mono or stereo audio
%                (sound pressure) signals
%
% sampleRateTarget : integer
%                    the sample rate (frequency) of the input target signal(s)
%
% axisTarget : integer (1 or 2, default: 1)
%              the time axis for the target signal(s) along which to determine
%              detection
%
% signalMasker : vector or 2D matrix
%                the input masker signal(s) as single mono or stereo audio
%                (sound pressure) signals
%
% sampleRateMasker : integer
%                    the sample rate (frequency) of the input masker signal(s)
%
% axisMasker : integer (1 or 2, default: 1)
%              the time axis for the masker signal(s) along which to determine
%              detection
%
% timeStep : number (default: 0.5)
%            the time window (seconds) to use for calculating target
%            detectability
%
% freqBandRange : vector (default: [20, 20000])
%                 the 1/3-octave band range over which to determine
%                 detection and discounted spectra
%
% outPlot : Boolean (default: false)
%           determines whether to plot outputs from the calculations
% 
% Returns
% -------
% detectDiscount : structure
%                  contains the output
%
% detectDiscount contains the following outputs:
%
% 
%
% Assumptions
% -----------
% The input signals are calibrated to units of acoustic pressure in Pascals
% (Pa).
%
% Requirements
% ------------
% Signal Processing Toolbox
% Audio Toolbox
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 05/11/2024
% Date last modified: 05/11/2024
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
        signalTarget (:, :) double {mustBeReal}
        sampleRateTarget (1, 1) double {mustBePositive, mustBeInteger}
        signalMasker (:, :) double {mustBeReal}
        sampleRateMasker (1, 1) double {mustBePositive, mustBeInteger}
        axisTarget (1, 1) {mustBeInteger, mustBeInRange(axisTarget, 1, 2)} = 1
        axisMasker (1, 1) {mustBeInteger, mustBeInRange(axisMasker, 1, 2)} = 1
        timeStep (1, 1) double {mustBePositive} = 0.5
        freqBandRange (1, 2) double {mustBeInRange(freqBandRange, 16, 22000)} = [20, 20000]
        outPlot {mustBeNumericOrLogical} = false
    end

%% Input checks and resampling


%% Define constants
b = 3;  % denominator for 1/b-octave definition

timeSteps = timeStep*sampleRate;

efficiencyFactor= 0.3;  % \eta
targetDetect = 2;  % target d'
discountHalfPower = 3;  % \alpha, dB
detectKnee = 14;  % \delta, d' value at which 3 dB knee occurs
discountRate = 1;  % \rho, rate at which discount function dimishes with reducing detectability

fl = min(freqBandRange);
fh = max(freqBandRange);

% free-field frontal tone hearing thresholds for 18-25 year-olds with
% normal hearing from ISO 226:2023 (20 Hz - 12.5 kHz)
hearThresholds226 = [78.1; 68.7; 59.5; 51.1; 44.0; 37.5; 31.5; 26.5; 22.1;...
                     17.9; 14.4; 11.4; 8.6; 6.2; 4.4; 3.0; 2.2; 2.4; 3.5;...
                     1.7; -1.3; -4.2; -6.0; -5.4; -1.5; 6.0; 12.6; 13.9; 12.3];

% free-field frontal tone hearing thresholds for 18-25 year-olds with
% normal hearing from ISO 389-7:2019 (14 - 18 kHz)
hearThresholds3897 = [18.4; 40.2; 70.4];

% free-field frontal tone hearing thresholds for 18-25 year-olds with
% normal hearing from ISO 226:2023 | ISO 389-7:2019 (20 - 18 kHz), with 20
% kHz band estimated from figure 1 in ISO 389-7:2019
hearThresholds = [hearThresholds226; hearThresholds3897; 100];

% diffuse field narrowband noise hearing thresholds for 18-25 year-olds with
% normal hearing from ISO 389-7:2019 (20 Hz - 16 kHz)
% NOTE: watch out for odd frequencies in the standard Table 1!
hearThresholdsDF3897 = [78.1; 68.7; 59.5; 51.1; 44.0; 37.5; 31.5; 26.5; 22.1;...
                        17.9; 14.4; 11.4; 8.4; 5.8; 3.8; 2.1; 1.0; 0.8; 1.9;...
                        0.5; -1.5; -3.1; -4.0; -3.8; -1.8; 2.5; 6.8; 8.4; 14.4;...
                        23.2; 43.7];

% estimated full range diffuse field narrowband noise hearing thresholds
% for 18-25 year-olds with normal hearing 
hearThresholdsDF = [hearThresholdsDF3897; 70.4; 100];

% 1/3-octave A-weighting dB values (20 Hz - 20 kHz)
Aweight = [-50.5; -44.7; -39.4; -34.6; -30.2; -26.2; -22.5; -19.1; -16.1;...
           -13.4; -10.9; -8.6; -6.6; -4.8; -3.2; -1.9; -0.8; 0.0; 0.6; 1.0;...
            1.2; 1.3; 1.2; 1.0; 0.5; -0.1; -1.1; -2.5; -4.3; -6.6; -9.3];


%% Signal processing

% Calculate 1/3-octave power spectrograms

[spectroTarget, ~, ~] = poctave(signalTarget, sampleRate, 'spectrogram',...
                                'BandsPerOctave', b, 'WindowLength', timeSteps,...
                                'FrequencyLimits', [fl, fh]);

[spectroMasker, f, t] = poctave(signalMasker, sampleRate, 'spectrogram',...
                                'BandsPerOctave', b, 'WindowLength', timeSteps,...
                                'FrequencyLimits', [fl, fh]);

% Equivalent rectangular bandwidth for auditory filter (Glasberg & Moore,
% 1990, equation 3)
ERBandwidth = 24.7*(1 + 4.37/1000*f);

% Calculate 1/3 octave bandwidths according to IEC 61640-1:2014 / ANSI S1.11-2004
G10 = 10^(3/10);  % octave ratio coefficient (base-ten)
OctRatio = G10^(0.5/b);  % octave ratio

f1 = f/OctRatio;  % output range of exact lower band-edge frequencies
f2 = f*OctRatio;  % output range of exact upper band-edge frequencies
fBandWidth = f2 - f1;

% Calculate high frequency weighting "kick" factor
wKick = zeros(size(f));
wKick(f >= 2500) = -6*(log10(f( f>= 2500)/2500)).^2;

% Calculate detection efficiency
detectEfficiency = efficiencyFactor*sqrt(timeStep)*sqrt(fBandWidth).*sqrt(fBandWidth./ERBandwidth).*10.^(wKick/10);

% Calculate equivalent auditory system noise
ind = (-17:1:13).';  % range of frequency indices for hearing threshold bands
if mod(b, 1) == 0
    fm = G10.^(ind/b)*1000;
else
    fm = G10.^((2*ind + 1)/(2*b))*1000;
end
[~, il] = min(abs(fm - f(1)));  % find nearest lower exact frequency
[~, ih] = min(abs(fm - f(end)));  % find nearest higher exact frequency

eqAuditoryNoise = detectEfficiency.*(2e-5.*10.^(hearThresholdsDF(il:ih)/20)).^2/targetDetect;

% Calculate detectability and discount
detectability = detectEfficiency.*spectroTarget./(spectroMasker + eqAuditoryNoise);
discountdB = discountHalfPower./(detectability./detectKnee).^discountRate;

% Calculate discounted levels
dBSpecTarget = 20*log10(sqrt(spectroTarget)/2e-5);
dBSpecMasker = 20*log10(sqrt(spectroMasker)/2e-5);
dBSpecDiscTarget = dBSpecTarget - discountdB;

% A-weight time-dependent spectra
dBASpecTarget = dBSpecTarget + Aweight;
dBASpecMasker = dBSpecMasker + Aweight;
dBASpecDiscTarget = dBSpecDiscTarget + Aweight;

% A-weighted time-dependent spectral power
powATarget = (2e-5*10.^(dBASpecTarget/20)).^2;
powAMasker = (2e-5*10.^(dBASpecMasker/20)).^2;
powADiscTarget = (2e-5*10.^(dBASpecDiscTarget/20)).^2;

% Aggregate spectral levels
dBATDepTarget = squeeze(20*log10(sqrt(sum(powATarget, 1))/2e-5));
dBATDepMasker = squeeze(20*log10(sqrt(sum(powAMasker, 1))/2e-5));
dBATDepDiscTarget = squeeze(20*log10(sqrt(sum(powADiscTarget, 1))/2e-5));

% Aggregate time-dependent levels
LAETarget = squeeze(20*log10(sqrt(sum(sum(powATarget, 1), 2))/2e-5));
LAEMasker = squeeze(20*log10(sqrt(sum(sum(powAMasker, 1), 2))/2e-5));
LAEDiscTarget = squeeze(20*log10(sqrt(sum(sum(powADiscTarget, 1), 2))/2e-5));

%% Assign outputs
detectDiscount.dBSpecTarget = dBSpecTarget;
detectDiscount.dBSpecMasker = dBSpecMasker;
detectDiscount.dBSpecDiscTarget = dBSpecDiscTarget;
detectDiscount.dBTDepTarget = dBATDepTarget;
detectDiscount.dBTDepMasker = dBATDepMasker;
detectDiscount.dBTDepDiscTarget = dBATDepDiscTarget;
detectDiscount.LAETarget = LAETarget;
detectDiscount.LAEMasker = LAEMasker;
detectDiscount.LAEDiscTarget = LAEDiscTarget;
detectDiscount.detectability = detectability;
detectDiscount.freqBands = f;
detectDiscount.timeOut = t;

%% Plotting

% end of function