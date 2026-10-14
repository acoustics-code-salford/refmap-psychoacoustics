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
addpath(genpath(fullfile(pwd)))

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

metricCount = 12;
% unmodulated high-pass noise
hprUnModSQMs = zeros(size(hprUnMod, 2), metricCount);
for ii = 1:size(hprUnMod, 2)
    PA = acousticPsychAnnoy(hprUnMod(:, ii), fs, 1, 2, "freeFrontal");
    hprUnModSQMs(ii, 1) = PA.loudness.N5;
    hprUnModSQMs(ii, 2) = PA.sharpness.S5;
    hprUnModSQMs(ii, 3) = PA.fluctuation.F5;
    hprUnModSQMs(ii, 4) = PA.roughness.R5;
    hprUnModSQMs(ii, 5) = PA.sharpness.Smean;
    hprUnModSQMs(ii, 6) = PA.fluctuation.Fmean;
    hprUnModSQMs(ii, 7) = PA.roughness.Rmean;
    hprUnModSQMs(ii, 8) = PA.fluctuation.F10;
    hprUnModSQMs(ii, 9) = PA.roughness.R10;
    hprUnModSQMs(ii, 10) = PA.psychAnnoyFrom5Pc;
    hprUnModSQMs(ii, 11) = prctile(PA.psychAnnoyTDep, 95);
    hprUnModSQMs(ii, 12) = PA.psychAnnoyFromMean;
end

% modulated high-pass noise
hprModSQMs = zeros(size(hprMod, 2), metricCount);
for ii = 1:size(hprMod, 2)
    PA = acousticPsychAnnoy(hprMod(:, ii), fs, 1, 2, "freeFrontal");
    hprModSQMs(ii, 1) = PA.loudness.N5;
    hprModSQMs(ii, 2) = PA.sharpness.S5;
    hprModSQMs(ii, 3) = PA.fluctuation.F5;
    hprModSQMs(ii, 4) = PA.roughness.R5;
    hprModSQMs(ii, 5) = PA.sharpness.Smean;
    hprModSQMs(ii, 6) = PA.fluctuation.Fmean;
    hprModSQMs(ii, 7) = PA.roughness.Rmean;
    hprModSQMs(ii, 8) = PA.fluctuation.F10;
    hprModSQMs(ii, 9) = PA.roughness.R10;
    hprModSQMs(ii, 10) = PA.psychAnnoyFrom5Pc;
    hprModSQMs(ii, 11) = prctile(PA.psychAnnoyTDep, 95);
    hprModSQMs(ii, 12) = PA.psychAnnoyFromMean;
end

% unmodulated blue noise
blrUnModSQMs = zeros(size(blrUnMod, 2), metricCount);
for ii = 1:size(blrUnMod, 2)
    PA = acousticPsychAnnoy(blrUnMod(:, ii), fs, 1, 2, "freeFrontal");
    blrUnModSQMs(ii, 1) = PA.loudness.N5;
    blrUnModSQMs(ii, 2) = PA.sharpness.S5;
    blrUnModSQMs(ii, 3) = PA.fluctuation.F5;
    blrUnModSQMs(ii, 4) = PA.roughness.R5;
    blrUnModSQMs(ii, 5) = PA.sharpness.Smean;
    blrUnModSQMs(ii, 6) = PA.fluctuation.Fmean;
    blrUnModSQMs(ii, 7) = PA.roughness.Rmean;
    blrUnModSQMs(ii, 8) = PA.fluctuation.F10;
    blrUnModSQMs(ii, 9) = PA.roughness.R10;
    blrUnModSQMs(ii, 10) = PA.psychAnnoyFrom5Pc;
    blrUnModSQMs(ii, 11) = prctile(PA.psychAnnoyTDep, 95);
    blrUnModSQMs(ii, 12) = PA.psychAnnoyFromMean;
end

% modulated blue noise
blrModSQMs = zeros(size(blrMod, 2), metricCount);
for ii = 1:size(blrMod, 2)
    PA = acousticPsychAnnoy(blrMod(:, ii), fs, 1, 2, "freeFrontal");
    blrModSQMs(ii, 1) = PA.loudness.N5;
    blrModSQMs(ii, 2) = PA.sharpness.S5;
    blrModSQMs(ii, 3) = PA.fluctuation.F5;
    blrModSQMs(ii, 4) = PA.roughness.R5;
    blrModSQMs(ii, 5) = PA.sharpness.Smean;
    blrModSQMs(ii, 6) = PA.fluctuation.Fmean;
    blrModSQMs(ii, 7) = PA.roughness.Rmean;
    blrModSQMs(ii, 8) = PA.fluctuation.F10;
    blrModSQMs(ii, 9) = PA.roughness.R10;
    blrModSQMs(ii, 10) = PA.psychAnnoyFrom5Pc;
    blrModSQMs(ii, 11) = prctile(PA.psychAnnoyTDep, 95);
    blrModSQMs(ii, 12) = PA.psychAnnoyFromMean;
end


% unmodulated uniform masking noise
gvrUnModSQMs = zeros(size(gvrUnMod, 2), metricCount);
for ii = 1:size(gvrUnMod, 2)
    PA = acousticPsychAnnoy(gvrUnMod(:, ii), fs, 1, 2, "freeFrontal");
    gvrUnModSQMs(ii, 1) = PA.loudness.N5;
    gvrUnModSQMs(ii, 2) = PA.sharpness.S5;
    gvrUnModSQMs(ii, 3) = PA.fluctuation.F5;
    gvrUnModSQMs(ii, 4) = PA.roughness.R5;
    gvrUnModSQMs(ii, 5) = PA.sharpness.Smean;
    gvrUnModSQMs(ii, 6) = PA.fluctuation.Fmean;
    gvrUnModSQMs(ii, 7) = PA.roughness.Rmean;
    gvrUnModSQMs(ii, 8) = PA.fluctuation.F10;
    gvrUnModSQMs(ii, 9) = PA.roughness.R10;
    gvrUnModSQMs(ii, 10) = PA.psychAnnoyFrom5Pc;
    gvrUnModSQMs(ii, 11) = prctile(PA.psychAnnoyTDep, 95);
    gvrUnModSQMs(ii, 12) = PA.psychAnnoyFromMean;
end

% modulated uniform masking noise
gvrModSQMs = zeros(size(gvrMod, 2), metricCount);
for ii = 1:size(gvrMod, 2)
    PA = acousticPsychAnnoy(gvrMod(:, ii), fs, 1, 2, "freeFrontal");
    gvrModSQMs(ii, 1) = PA.loudness.N5;
    gvrModSQMs(ii, 2) = PA.sharpness.S5;
    gvrModSQMs(ii, 3) = PA.fluctuation.F5;
    gvrModSQMs(ii, 4) = PA.roughness.R5;
    gvrModSQMs(ii, 5) = PA.sharpness.Smean;
    gvrModSQMs(ii, 6) = PA.fluctuation.Fmean;
    gvrModSQMs(ii, 7) = PA.roughness.Rmean;
    gvrModSQMs(ii, 8) = PA.fluctuation.F10;
    gvrModSQMs(ii, 9) = PA.roughness.R10;
    gvrModSQMs(ii, 10) = PA.psychAnnoyFrom5Pc;
    gvrModSQMs(ii, 11) = prctile(PA.psychAnnoyTDep, 95);
    gvrModSQMs(ii, 12) = PA.psychAnnoyFromMean;
end

% calculate relative values
hprUnModSQMs(:, metricCount + 1) = hprUnModSQMs(:, 1)./gvrUnModSQMs(2, 1);
hprUnModSQMs(:, metricCount + 2) = hprUnModSQMs(:, metricCount - 2)./gvrUnModSQMs(2, metricCount - 2);
hprModSQMs(:, metricCount + 1) = hprModSQMs(:, 1)./gvrUnModSQMs(2, 1);
hprModSQMs(:, metricCount + 2) = hprModSQMs(:, metricCount - 2)./gvrUnModSQMs(2, metricCount - 2);
blrUnModSQMs(:, metricCount + 1) = blrUnModSQMs(:, 1)./gvrUnModSQMs(2, 1);
blrUnModSQMs(:, metricCount + 2) = blrUnModSQMs(:, metricCount - 2)./gvrUnModSQMs(2, metricCount - 2);
blrModSQMs(:, metricCount + 1) = blrModSQMs(:, 1)./gvrUnModSQMs(2, 1);
blrModSQMs(:, metricCount + 2) = blrModSQMs(:, metricCount - 2)./gvrUnModSQMs(2, metricCount - 2);
gvrUnModSQMs(:, metricCount + 1) = gvrUnModSQMs(:, 1)./gvrUnModSQMs(2, 1);
gvrUnModSQMs(:, metricCount + 2) = gvrUnModSQMs(:, metricCount - 2)./gvrUnModSQMs(2, metricCount - 2);
gvrModSQMs(:, metricCount + 1) = gvrModSQMs(:, 1)./gvrUnModSQMs(2, 1);
gvrModSQMs(:, metricCount + 2) = gvrModSQMs(:, metricCount - 2)./gvrUnModSQMs(2, metricCount - 2);


%% Print table
tab4_2 = vertcat(gvrUnModSQMs(1, :),...
    gvrModSQMs(1, :),...
    blrUnModSQMs(1, :),...
    blrModSQMs(1, :),...
    hprUnModSQMs(1, :),...
    hprModSQMs(1, :),...
    gvrUnModSQMs(2, :),...
    gvrModSQMs(2, :),...
    blrUnModSQMs(2, :),...
    blrModSQMs(2, :),...
    hprUnModSQMs(2, :),...
    hprModSQMs(2, :),...
    gvrUnModSQMs(3, :),...
    gvrModSQMs(3, :),...
    blrUnModSQMs(3, :),...
    blrModSQMs(3, :),...
    hprUnModSQMs(3, :),...
    hprModSQMs(3, :));


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