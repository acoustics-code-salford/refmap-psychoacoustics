function acousticSHMGenerateRefSignals
% acousticSHMGenerateRefSignals
%
% Generates reference signals as wav format mono audio files for
% calibrating and testing the validity of ECMA-418-2 implementation (Sottek
% Hearing Model)
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 19/08/2024
% Date last modified: 21/07/2025
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

%% Load path (assumes root directory is refmap-psychoacoustics)
addpath(genpath(fullfile("validation")))

%% Input parameters

fs = 48e3;
dt = 1/fs;
T = 10;
n = T*fs;
t = linspace(0, T-dt, n);
f_tone = 1000;
f_mod70 = 70;
f_mod4 = 4;
A_tone = sqrt(2)*2e-5*10^(40/20);
outpath = fullfile("validation",...
                   "ECMA-418-2", "audio");

%% Generate signals
% reference sinusoid for loudness and tonality

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