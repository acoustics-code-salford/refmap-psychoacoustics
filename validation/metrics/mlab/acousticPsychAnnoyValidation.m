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

% unmodulated high-pass noise
hprUnModSQMs = zeros(7, size(hprUnMod, 2));
for ii = 1:size(hprUnMod, 2)
    PA = acousticPsychAnnoy(hprUnMod(:, ii), fs, 1, 2, "freeFrontal");
    hprUnModSQMs(1, ii) = PA.loudness.N5;
    hprUnModSQMs(2, ii) = PA.sharpness.S5;
    hprUnModSQMs(3, ii) = PA.fluctuation.F5;
    hprUnModSQMs(4, ii) = PA.roughness.R5;
    hprUnModSQMs(5, ii) = PA.psychAnnoyFrom5Pc;
end

% modulated high-pass noise
hprModSQMs = zeros(7, size(hprMod, 2));
for ii = 1:size(hprMod, 2)
    PA = acousticPsychAnnoy(hprMod(:, ii), fs, 1, 2, "freeFrontal");
    hprModSQMs(1, ii) = PA.loudness.N5;
    hprModSQMs(2, ii) = PA.sharpness.S5;
    hprModSQMs(3, ii) = PA.fluctuation.F5;
    hprModSQMs(4, ii) = PA.roughness.R5;
    hprModSQMs(5, ii) = PA.psychAnnoyFrom5Pc;
end

% unmodulated blue noise
blrUnModSQMs = zeros(7, size(blrUnMod, 2));
for ii = 1:size(blrUnMod, 2)
    PA = acousticPsychAnnoy(blrUnMod(:, ii), fs, 1, 2, "freeFrontal");
    blrUnModSQMs(1, ii) = PA.loudness.N5;
    blrUnModSQMs(2, ii) = PA.sharpness.S5;
    blrUnModSQMs(3, ii) = PA.fluctuation.F5;
    blrUnModSQMs(4, ii) = PA.roughness.R5;
    blrUnModSQMs(5, ii) = PA.psychAnnoyFrom5Pc;
end

% modulated blue noise
blrModSQMs = zeros(7, size(blrMod, 2));
for ii = 1:size(blrMod, 2)
    PA = acousticPsychAnnoy(blrMod(:, ii), fs, 1, 2, "freeFrontal");
    blrModSQMs(1, ii) = PA.loudness.N5;
    blrModSQMs(2, ii) = PA.sharpness.S5;
    blrModSQMs(3, ii) = PA.fluctuation.F5;
    blrModSQMs(4, ii) = PA.roughness.R5;
    blrModSQMs(5, ii) = PA.psychAnnoyFrom5Pc;
end


% unmodulated uniform masking noise
gvrUnModSQMs = zeros(7, size(gvrUnMod, 2));
for ii = 1:size(gvrUnMod, 2)
    PA = acousticPsychAnnoy(gvrUnMod(:, ii), fs, 1, 2, "freeFrontal");
    gvrUnModSQMs(1, ii) = PA.loudness.N5;
    gvrUnModSQMs(2, ii) = PA.sharpness.S5;
    gvrUnModSQMs(3, ii) = PA.fluctuation.F5;
    gvrUnModSQMs(4, ii) = PA.roughness.R5;
    gvrUnModSQMs(5, ii) = PA.psychAnnoyFrom5Pc;
end

% modulated uniform masking noise
gvrModSQMs = zeros(7, size(gvrMod, 2));
for ii = 1:size(gvrMod, 2)
    PA = acousticPsychAnnoy(gvrMod(:, ii), fs, 1, 2, "freeFrontal");
    gvrModSQMs(1, ii) = PA.loudness.N5;
    gvrModSQMs(2, ii) = PA.sharpness.S5;
    gvrModSQMs(3, ii) = PA.fluctuation.F5;
    gvrModSQMs(4, ii) = PA.roughness.R5;
    gvrModSQMs(5, ii) = PA.psychAnnoyFrom5Pc;
end

% calculate relative values
hprUnModSQMs(6, :) = hprUnModSQMs(1, :)./gvrUnModSQMs(1, 2);
hprUnModSQMs(7, :) = hprUnModSQMs(5, :)./gvrUnModSQMs(5, 2);
hprModSQMs(6, :) = hprModSQMs(1, :)./gvrUnModSQMs(1, 2);
hprModSQMs(7, :) = hprModSQMs(5, :)./gvrUnModSQMs(5, 2);
blrUnModSQMs(6, :) = blrUnModSQMs(1, :)./gvrUnModSQMs(1, 2);
blrUnModSQMs(7, :) = blrUnModSQMs(5, :)./gvrUnModSQMs(5, 2);
blrModSQMs(6, :) = blrModSQMs(1, :)./gvrUnModSQMs(1, 2);
blrModSQMs(7, :) = blrModSQMs(5, :)./gvrUnModSQMs(5, 2);
gvrUnModSQMs(6, :) = gvrUnModSQMs(1, :)./gvrUnModSQMs(1, 2);
gvrUnModSQMs(7, :) = gvrUnModSQMs(5, :)./gvrUnModSQMs(5, 2);
gvrModSQMs(6, :) = gvrModSQMs(1, :)./gvrUnModSQMs(1, 2);
gvrModSQMs(7, :) = gvrModSQMs(5, :)./gvrUnModSQMs(5, 2);

%% Plot results

% arrange data
relLoudData1 = [gvrUnModSQMs(6, 1), gvrModSQMs(6, 1),...
                blrUnModSQMs(6, 1), blrModSQMs(6, 1),...
                hprUnModSQMs(6, 1), hprModSQMs(6, 1)];
relLoudData2 = [gvrUnModSQMs(6, 2), gvrModSQMs(6, 2),...
                blrUnModSQMs(6, 2), blrModSQMs(6, 2),...
                hprUnModSQMs(6, 2), hprModSQMs(6, 2)];
relLoudData3 = [gvrUnModSQMs(6, 3), gvrModSQMs(6, 3),...
                blrUnModSQMs(6, 3), blrModSQMs(6, 3),...
                hprUnModSQMs(6, 3), hprModSQMs(6, 3)];
relPAData1 = [gvrUnModSQMs(7, 1), gvrModSQMs(7, 1),...
              blrUnModSQMs(7, 1), blrModSQMs(7, 1),...
              hprUnModSQMs(7, 1), hprModSQMs(7, 1)];
relPAData2 = [gvrUnModSQMs(7, 2), gvrModSQMs(7, 2),...
              blrUnModSQMs(7, 2), blrModSQMs(7, 2),...
              hprUnModSQMs(7, 2), hprModSQMs(7, 2)];
relPAData3 = [gvrUnModSQMs(7, 3), gvrModSQMs(7, 3),...
              blrUnModSQMs(7, 3), blrModSQMs(7, 3),...
              hprUnModSQMs(7, 3), hprModSQMs(7, 3)];
plotLabels = {"uGVR", "mGVR", "uBLR", "mBLR", "uGVR", "mGVR"};

% relative loudness
fig = figure;
set(fig, 'position', [300, 200, 800, 350])
tiledlayout(fig, 1, 3);
movegui(fig, 'center');

ax1 = nexttile(1);
plot(1:6, relLoudData1, 'LineStyle', 'none', 'Marker', '^')
ax1.YLim = [0, 3];
ax1.XTick = [1:6];
ax1.XTickLabels = plotLabels;
ax1.XLim = [0, 7];
ax1.YLabel.String = "Relative loudness";

ax1 = nexttile(2);
plot(1:6, relLoudData2, 'LineStyle', 'none', 'Marker', '^')
ax1.YLim = [0, 3];
ax1.XTick = [1:6];
ax1.XTickLabels = plotLabels;
ax1.XLim = [0, 7];

ax1 = nexttile(3);
plot(1:6, relLoudData3, 'LineStyle', 'none', 'Marker', '^')
ax1.YLim = [0, 3];
ax1.XTick = [1:6];
ax1.XTickLabels = plotLabels;
ax1.XLim = [0, 7];


% relative psychoacoustic annoyance
fig = figure;
set(fig, 'position', [300, 200, 800, 500])
tiledlayout(fig, 1, 3);
movegui(fig, 'center');

ax1 = nexttile(1);
plot(1:6, relPAData1, 'LineStyle', 'none', 'Marker', 'x')
ax1.YLim = [0, 6];
ax1.XTick = [1:6];
ax1.XTickLabels = plotLabels;
ax1.XLim = [0, 7];
ax1.YLabel.String = "Relative psychoacoustic annoyance";

ax1 = nexttile(2);
plot(1:6, relPAData2, 'LineStyle', 'none', 'Marker', 'x')
ax1.YLim = [0, 6];
ax1.XTick = [1:6];
ax1.XTickLabels = plotLabels;
ax1.XLim = [0, 7];

ax1 = nexttile(3);
plot(1:6, relPAData3, 'LineStyle', 'none', 'Marker', 'x')
ax1.YLim = [0, 6];
ax1.XTick = [1:6];
ax1.XTickLabels = plotLabels;
ax1.XLim = [0, 7];


% end of function