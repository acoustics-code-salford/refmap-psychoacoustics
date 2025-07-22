function acousticPsychAnnoyGenerateRefSignals
% acousticPsychAnnoyGenerateRefSignals
%
% Generates reference signals as wav format mono audio files for
% calibrating and testing the validity of Widmann Psychoacoustic Annoyance
% (also misattributed as Zwicker PA). The 'reference signal' for the model
% was defined as a 1 kHz sinusoid at 40 dB Lp in free-field conditions. The
% signals generated here correspond with Table 4.2 of Widmann (1992), which
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
% Date last modified: 13/07/2025
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

fs = 48e3;  % sampling frequency
dt = 1/fs;  % time step
T = 10;  % total time
N = T*fs;  % total sample length
t = linspace(0, T - dt, N).';  % time vector
f_mod32 = 32;  % modulation frequency
m_mod = 0.98;  % modulation factor
outpath = fullfile("validation",...
                   "metrics", "audio");

%% Generate signals

% signal parameters are assigned according to Section 4.3 Table 4.2
% high-pass white noise (HPR)
whiteNoise = dsp.ColoredNoise('white', N, 1);
whiteNoiseSig = bandpass(whiteNoise(), [4e3, 20e3], fs);

% blue noise
blueNoise = dsp.ColoredNoise('blue', N, 1);
blueNoiseSig = bandpass(blueNoise(), [20, 20e3], fs);

% uniform masking noise
[unifMaskNoise, ~] = audioread("uniform_masking_noise.wav");

% modulation signal

sine_mod_32Hz_40dB = m_mod*sin(2*pi*f_mod32.*t);

%% Generate signals and save as wav files
% (assumes current folder is refmap-psychoacoustics root)

% high-pass noise
% unmodulated signals
whiteUnModLoud = [4.4, 10.5, 22.6];
hprUnMod = zeros(N, length(whiteUnModLoud));
for ii = 1:length(whiteUnModLoud)
    
    hprUnMod(:, ii) = setLoudness(whiteNoiseSig, fs, whiteUnModLoud(ii),...
                                  1, "freeFrontal", "5PcEx",...
                                  "ISO 532-1", false, 0.01, 15);
    filename = "hpassWNoise_4kHz_N" + num2str(round(whiteUnModLoud(ii))) + ".wav";

    % select required bit depth
    if max(abs(hprUnMod(:, ii))) > 1
        bitsRequired = 32;
    else
        bitsRequired = 24;
    end

    audiowrite(fullfile(outpath, filename), hprUnMod(:, ii), fs,...
               "BitsPerSample", bitsRequired)
end

% modulated signals
whiteNoiseSigMod = (1 + sine_mod_32Hz_40dB).*whiteNoiseSig;
whiteModLoud = [4.4, 10.6, 22.8];
hprMod = zeros(N, length(whiteModLoud));
for ii = 1:length(whiteModLoud)
    hprMod(:, ii) = setLoudness(whiteNoiseSigMod, fs, whiteModLoud(ii),...
                                1, "freeFrontal", "5PcEx",...
                                "ISO 532-1", false, 0.01, 15);
    filename = "hpassWNoise_4kHz_N" + num2str(round(whiteModLoud(ii))) + "_mod_" + num2str(f_mod32) + "Hz.wav";

    % select required bit depth
    if max(abs(hprUnMod(:, ii))) > 1
        bitsRequired = 32;
    else
        bitsRequired = 24;
    end

    audiowrite(fullfile(outpath, filename), hprMod(:, ii), fs,...
               "BitsPerSample", bitsRequired)
end

% blue noise
% unmodulated signals
blueUnModLoud = [4.2, 9.4, 19.9];
blrUnMod = zeros(N, length(blueUnModLoud));
for ii = 1:length(blueUnModLoud)
    blrUnMod(:, ii) = setLoudness(blueNoiseSig, fs, blueUnModLoud(ii),...
                                  1, "freeFrontal", "5PcEx",...
                                  "ISO 532-1", false, 0.01, 15);
    filename = "blueNoise_N" + num2str(round(blueUnModLoud(ii))) + ".wav";

    % select required bit depth
    if max(abs(hprUnMod(:, ii))) > 1
        bitsRequired = 32;
    else
        bitsRequired = 24;
    end

    audiowrite(fullfile(outpath, filename), blrUnMod(:, ii), fs,...
               "BitsPerSample", bitsRequired)
end

% modulated signals
blueNoiseSigMod = (1 + sine_mod_32Hz_40dB).*blueNoiseSig;
blueModLoud = [4.3, 9.5, 20.5];
blrMod = zeros(N, length(blueModLoud));
for ii = 1:length(blueModLoud)
    blrMod(:, ii) = setLoudness(blueNoiseSigMod, fs, blueModLoud(ii),...
                                1, "freeFrontal", "5PcEx",...
                                "ISO 532-1", false, 0.01, 15);
    filename = "blueNoise_N" + num2str(round(blueModLoud(ii))) + "_mod_" + num2str(f_mod32) + "Hz.wav";

    % select required bit depth
    if max(abs(hprUnMod(:, ii))) > 1
        bitsRequired = 32;
    else
        bitsRequired = 24;
    end

    audiowrite(fullfile(outpath, filename), blrMod(:, ii), fs,...
               "BitsPerSample", bitsRequired)
end

% uniform masking noise
% unmodulated signals
unifMaskUnModLoud = [4.5, 10.8, 23];
gvrUnMod = zeros(N, length(unifMaskUnModLoud));
for ii = 1:length(unifMaskUnModLoud)
    gvrUnMod(:, ii) = setLoudness(unifMaskNoise, fs, unifMaskUnModLoud(ii),...
                                  1, "freeFrontal", "5PcEx",...
                                  "ISO 532-1", false, 0.01, 15);
    filename = "unifMaskNoise_N" + num2str(round(unifMaskUnModLoud(ii))) + ".wav";

    % select required bit depth
    if max(abs(hprUnMod(:, ii))) > 1
        bitsRequired = 32;
    else
        bitsRequired = 24;
    end

    audiowrite(fullfile(outpath, filename), gvrUnMod(:, ii), fs,...
               "BitsPerSample", bitsRequired)
end

% modulated signals
unifMaskNoiseSigMod = (1 + sine_mod_32Hz_40dB).*unifMaskNoise;
unifMaskModLoud = [4.6, 10.9, 22.8];
gvrMod = zeros(N, length(unifMaskModLoud));
for ii = 1:length(unifMaskModLoud)
    gvrMod(:, ii) = setLoudness(unifMaskNoiseSigMod, fs, unifMaskModLoud(ii),...
                                1, "freeFrontal", "5PcEx",...
                                "ISO 532-1", false, 0.01, 15);
    filename = "unifMaskNoise_N" + num2str(round(unifMaskModLoud(ii))) + "_mod_" + num2str(f_mod32) + "Hz.wav";

    % select required bit depth
    if max(abs(hprUnMod(:, ii))) > 1
        bitsRequired = 32;
    else
        bitsRequired = 24;
    end

    audiowrite(fullfile(outpath, filename), gvrMod(:, ii), fs,...
               "BitsPerSample", bitsRequired)
end

% end of function