function acousticHMSGenerateRefSignals
% acousticHMSGenerateRefSignals
%
% Generates reference signals as wav format mono audio files for
% calibrating and testing the validity of ECMA-418-2:2022 implementation.
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 19/08/2024
% Date last modified: 19/08/2024
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
% input parameters

fs = 48e3;
dt = 1/fs;
T = 5;
n = T*fs;
t = linspace(0, T-dt, n);
f_tone = 1000;
f_mod = 70;
A_tone = sqrt(2)*2e-5*10^(40/20);
outpath = fullfile("validation",...
                   "ECMA-418-2_2022", "audio");

% reference sinusoid for loudness and tonality

sine_1kHz_40dB = A_tone*sin(2*pi*f_tone.*t);

% reference modulated sinusoid for roughness

sine_70Hz_mod = sin(2*pi*f_mod.*t - pi/2);
sine_1kHz_70Hz_60dB = (1 + sine_70Hz_mod).*sin(2*pi*f_tone.*t);
A_adjust = rms(sine_1kHz_70Hz_60dB)/0.02;
sine_1kHz_70Hz_60dB = sine_1kHz_70Hz_60dB/A_adjust;

% save signals as wav files (assumes current folder is
% refmap-psychoacoustics root)
audiowrite(fullfile(outpath, "sine_1kHz_40dB.wav"), sine_1kHz_40dB, fs,...
           "BitsPerSample", 24)
audiowrite(fullfile(outpath, "sine_1kHz_70Hz_60dB.wav"), sine_1kHz_70Hz_60dB, fs,...
           "BitsPerSample", 24)

% end of function