function acousticPsychAnnoyGenerateRefSignals
% acousticPsychAnnoyGenerateRefSignals
%
% Generates reference signals as wav format mono audio files for
% calibrating and testing the validity of Widmann Psychoacoustic Annoyance
% (also misattributed as Zwicker PA). The 'reference signal' for the model
% was defined as a 1 kHz sinusoid at 40 dB Lp in free-field conditions. The
% signals generated here correspond with Table 4.1 of Widmann (1992), which
% offer a wider range of sound characteristics to compare with calculated
% PA values.
%
% Widmann, U. (1992). Ein Modell der Psychoakustischen Lästigkeit von
% Schallen und seine Anwendung in der Praxis der Lärmbeurteilung (A model
% of the psychoacoustic annoyance of sounds and its application in noise
% assessment practice) [Doctoral thesis, Technische Universität München
% (Technical University of Munich)].
%
% Commonly misattributed to Zwicker, E. and Fastl, H. Second ed.,
% Psychoacoustics, Facts and Models, 2nd ed. Springer-Verlag, Berlin, 1999.
%
% Lotinga, M. J. B. and A. J. Torija (2025). "Comment on "A study on
% calibration methods of noise annoyance data from listening tests"
% [J. Acoust. Soc. Am. 156, 1877–1886 (2024)]." Journal of the Acoustical
% Society of America 157(5): 3282–3285.
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 10/07/2025
% Date last modified: 11/07/2025
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

%% Load path
addpath(genpath(fullfile("refmap-psychoacoustics")))

%% Input parameters

seed = 808;
fs = 48e3;
dt = 1/fs;
T = 10;
N = T*fs;
t = linspace(0, T - dt, N);
f_mod32 = 32;
f_mod1 = 1;
outpath = fullfile("validation",...
                   "metrics", "audio");

%% Generate signals

% signal parameters are assigned according to Table 4.1
% white noise-based
whiteNoise = dsp.ColoredNoise('white', N, 1);
whiteNoiseSig = bandpass(whiteNoise(), [20, 20e3], fs);

% pink noise-based
pinkNoise = dsp.ColoredNoise('pink', N, 1);
pinkNoiseSig = bandpass(pinkNoise(), [20, 20e3], fs);

% blue noise-based
blueNoise = dsp.ColoredNoise('blue', N, 1);
blueNoiseSig = bandpass(blueNoise(), [20, 20e3], fs);

sine_1kHz_40dB = A_tone*sin(2*pi*f_tone.*t);

% reference modulated sinusoid for roughness

sine_70Hz_mod = sin(2*pi*f_mod70.*t - pi/2);
sine_1kHz_70Hz_60dB = (1 + sine_70Hz_mod).*sin(2*pi*f_tone.*t);
A_adjust = rms(sine_1kHz_70Hz_60dB)/0.02;
sine_1kHz_70Hz_60dB = sine_1kHz_70Hz_60dB/A_adjust;

% reference modulated sinusoid for fluctuation strength

sine_4Hz_mod = sin(2*pi*f_mod4.*t - pi/2);
sine_1kHz_4Hz_60dB = (1 + sine_4Hz_mod).*sin(2*pi*f_tone.*t);
A_adjust = rms(sine_1kHz_4Hz_60dB)/0.02;
sine_1kHz_4Hz_60dB = sine_1kHz_4Hz_60dB/A_adjust;

%% Save signals as wav files (assumes current folder is
% refmap-psychoacoustics root)
audiowrite(fullfile(outpath, "sine_1kHz_40dB.wav"), sine_1kHz_40dB, fs,...
           "BitsPerSample", 24)
audiowrite(fullfile(outpath, "sine_1kHz_70Hz_60dB.wav"), sine_1kHz_70Hz_60dB, fs,...
           "BitsPerSample", 24)
audiowrite(fullfile(outpath, "sine_1kHz_4Hz_60dB.wav"), sine_1kHz_4Hz_60dB, fs,...
           "BitsPerSample", 24)

% end of function