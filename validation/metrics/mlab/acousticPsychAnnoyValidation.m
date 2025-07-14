function acousticPsychAnnoyValidation(savePlots)
% acousticPsychAnnoyValidation(savePlots)
%
% Compares Widmann psychoacoustic annoyance model implementation with 
% calculated values from Widmann (1992) for test signals defined therein.
% (Widmann PA is commonly misattributed as Zwicker PA).
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
% Inputs
% ------
% savePlots : Boolean (default = true)
%             flag indicating whether to save the generated plots to file
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 13/07/2025
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
%% Arguments validation
    arguments (Input)
        savePlots {mustBeNumericOrLogical} = true
    end

%% Load path
addpath(genpath(fullfile("refmap-psychoacoustics")))

%% Load test signals

% high-pass noise
% unmodulated signals
[hprUnMod(:, 1), fs] = audioread("hpassWNoise_4kHz_N4.wav");
[hprUnMod(:, 2), ~] = audioread("hpassWNoise_4kHz_N11.wav");
[hprUnMod(:, 3), ~] = audioread("hpassWNoise_4kHz_N23.wav");

% modulated signals
[hprMod(:, 1), ~] = audioread("hpassWNoise_4kHz_N4_mod_32Hz.wav");
[hprMod(:, 2), ~] = audioread("hpassWNoise_4kHz_N11_mod_32Hz.wav");
[hprMod(:, 3), ~] = audioread("hpassWNoise_4kHz_N23_mod_32Hz.wav");

% blue noise
% unmodulated signals
[blrUnMod(:, 1), ~] = audioread("blueNoise_N4.wav");
[blrUnMod(:, 2), ~] = audioread("blueNoise_N9.wav");
[blrUnMod(:, 3), ~] = audioread("blueNoise_N20.wav");

% modulated signals
[blrMod(:, 1), ~] = audioread("blueNoise_N4_mod_32Hz.wav");
[blrMod(:, 2), ~] = audioread("blueNoise_N10_mod_32Hz.wav");
[blrMod(:, 3), ~] = audioread("blueNoise_N21_mod_32Hz.wav");

% uniform masking noise
% unmodulated signals
[gvrUnMod(:, 1), ~] = audioread("unifMaskNoise_N5.wav");
[gvrUnMod(:, 2), ~] = audioread("unifMaskNoise_N11.wav");
[gvrUnMod(:, 3), ~] = audioread("unifMaskNoise_N23.wav");

% modulated signals
[gvrMod(:, 1), ~] = audioread("unifMaskNoise_N5_mod_32Hz.wav");
[gvrMod(:, 2), ~] = audioread("unifMaskNoise_N11_mod_32Hz.wav");
[gvrMod(:, 3), ~] = audioread("unifMaskNoise_N23_mod_32Hz.wav");

%% Comparison data from Widmann (1992)



%% Calculate predicted values for test signals

% end of function